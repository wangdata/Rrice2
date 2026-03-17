#!/usr/bin/env Rscript

# 再生稻头季-再生季产量关系综合分析脚本
# 功能：
# 1) 读取并清洗原始Excel数据
# 2) 基于经纬度补充土壤属性（优先 SoilGrids REST，可替换为本地 HWSD 栅格提取）
# 3) 补充头季/再生季气象变量（NASA POWER）
# 4) 连续变量分级
# 5) 缺失值、分布与相关性探索
# 6) 构建稳健模型（LMM + 鲁棒回归 + LASSO + 随机森林）

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(httr2)
  library(jsonlite)
  library(lubridate)
  library(nasapower)
  library(lme4)
  library(lmerTest)
  library(MASS)
  library(glmnet)
  library(ranger)
  library(performance)
  library(broom.mixed)
  library(car)
})

options(stringsAsFactors = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)

#-----------------------------
# 0. 参数区
#-----------------------------
input_file <- "20260201再生稻产量构成数据Meta.xlsx"
soil_cache <- "results/soil_cache.rds"
climate_cache <- "results/climate_cache.rds"

#-----------------------------
# 1. 读取与清洗
#-----------------------------
raw <- read_excel(input_file) |>
  clean_names() |>
  rename(
    author_id = author_id,
    data_id = data_id,
    year = year,
    lat = lati,
    lon = long,
    altitude = alti,
    country = country,
    variety = variety,
    genotype = genotype,
    cut_ht = cut_ht,
    n_main_crop = n_main_crop,
    n_ratoon = n_ratoon,
    n_total = n_total,
    m_panicle = mpanicle,
    m_spikelet = mspikelet,
    m_grain_wt = mgrain_wt,
    m_filled_rate = mfilled_grain_rate,
    m_yield = mgrain_yield,
    m_growth_length = mgrowth_length,
    r_panicle = rpanicle,
    r_spikelet = rspikelet,
    r_grain_wt = rgrain_wt,
    r_filled_rate = rfilled_grain_rate,
    r_yield = rgrain_yield,
    r_growth_length = rgrowth_length
  ) |>
  mutate(
    across(c(year, lat, lon, altitude, cut_ht, n_main_crop, n_ratoon, n_total,
             m_panicle, m_spikelet, m_grain_wt, m_filled_rate, m_yield, m_growth_length,
             r_panicle, r_spikelet, r_grain_wt, r_filled_rate, r_yield, r_growth_length),
           as.numeric),
    variety = as.factor(variety),
    genotype = as.factor(genotype),
    country = as.factor(country)
  )

#-----------------------------
# 2. 土壤属性补充（SoilGrids REST）
#-----------------------------
query_soilgrids <- function(lat, lon) {
  # ISRIC SoilGrids v2 REST API
  # 这里提取表层(0-5 cm)均值，可按研究需求替换深度层
  url <- paste0(
    "https://rest.isric.org/soilgrids/v2.0/properties/query",
    "?lon=", lon,
    "&lat=", lat,
    "&property=clay&property=nitrogen&property=soc&property=phh2o",
    "&depth=0-5cm",
    "&value=mean"
  )

  req <- request(url) |> req_user_agent("ratoon-rice-analysis/1.0")
  resp <- req_perform(req)
  txt <- resp_body_string(resp)
  js <- fromJSON(txt, flatten = TRUE)

  # 安全提取
  layers <- js$properties$layers
  get_val <- function(layer_name) {
    idx <- which(layers$name == layer_name)
    if (length(idx) == 0) return(NA_real_)
    as.numeric(layers$depths[[idx]][[1]]$values$mean)
  }

  clay <- get_val("clay")      # g/kg
  n_soil <- get_val("nitrogen")# cg/kg or dg/kg(数据源单位需按metadata核对)
  soc <- get_val("soc")        # dg/kg
  ph <- get_val("phh2o")       # pH*10

  tibble(
    lat = lat,
    lon = lon,
    soil_clay_gkg = clay,
    soil_n_raw = n_soil,
    soil_soc_raw = soc,
    soil_ph_raw = ph
  )
}

get_soil_data <- function(df) {
  locs <- df |> distinct(lat, lon) |> filter(!is.na(lat), !is.na(lon))

  if (file.exists(soil_cache)) {
    message("读取土壤缓存: ", soil_cache)
    soil <- readRDS(soil_cache)
    return(soil)
  }

  message("在线抓取土壤数据（SoilGrids API）...")
  soil <- purrr::pmap_dfr(locs, function(lat, lon) {
    tryCatch(query_soilgrids(lat, lon), error = function(e) {
      message("soil API失败 lat=", lat, ", lon=", lon, "; ", e$message)
      tibble(lat = lat, lon = lon,
             soil_clay_gkg = NA_real_, soil_n_raw = NA_real_,
             soil_soc_raw = NA_real_, soil_ph_raw = NA_real_)
    })
  })

  saveRDS(soil, soil_cache)
  soil
}

soil <- get_soil_data(raw) |>
  mutate(
    # SoilGrids变量常见缩放：SOC(dg/kg), pH*10（请以官方metadata复核）
    soil_soc_gkg = soil_soc_raw / 10,
    soil_som_gkg = soil_soc_gkg * 1.724, # SOM = SOC * 1.724
    soil_ph = soil_ph_raw / 10,
    soil_n_gkg = soil_n_raw / 100 # 示例换算，建议根据soilgrids单位再次确认
  )

#-----------------------------
# 3. 气象数据补充（头季/再生季分段）
#-----------------------------
infer_season_window <- function(year, lat, m_len, r_len) {
  # 无播期时的近似窗口：北半球默认4月1日开头季，南半球10月1日
  start_main <- ifelse(lat >= 0,
                       as.Date(sprintf("%d-04-01", year)),
                       as.Date(sprintf("%d-10-01", year)))
  end_main <- start_main + days(ifelse(is.na(m_len), 120, m_len) - 1)
  start_ratoon <- end_main + days(1)
  end_ratoon <- start_ratoon + days(ifelse(is.na(r_len), 60, r_len) - 1)
  tibble(start_main, end_main, start_ratoon, end_ratoon)
}

query_power_daily <- function(lat, lon, start_date, end_date) {
  # NASA POWER: T2M, PRECTOTCORR, ALLSKY_SFC_SW_DWN
  d <- tryCatch(
    get_power(
      community = "AG",
      lonlat = c(lon, lat),
      pars = c("T2M", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN"),
      dates = c(format(start_date, "%Y%m%d"), format(end_date, "%Y%m%d")),
      temporal_api = "DAILY"
    ),
    error = function(e) NULL
  )

  if (is.null(d) || nrow(d) == 0) {
    return(tibble(t2m = NA_real_, ppt = NA_real_, srad = NA_real_))
  }

  tibble(
    t2m = mean(d$T2M, na.rm = TRUE),
    ppt = sum(d$PRECTOTCORR, na.rm = TRUE),
    srad = mean(d$ALLSKY_SFC_SW_DWN, na.rm = TRUE)
  )
}

get_climate_data <- function(df) {
  key_df <- df |>
    distinct(data_id, year, lat, lon, m_growth_length, r_growth_length) |>
    filter(!is.na(year), !is.na(lat), !is.na(lon))

  if (file.exists(climate_cache)) {
    message("读取气象缓存: ", climate_cache)
    return(readRDS(climate_cache))
  }

  message("在线抓取气象数据（NASA POWER）...")
  out <- pmap_dfr(key_df, function(data_id, year, lat, lon, m_growth_length, r_growth_length) {
    w <- infer_season_window(year, lat, m_growth_length, r_growth_length)

    main_clim <- query_power_daily(lat, lon, w$start_main, w$end_main)
    ratoon_clim <- query_power_daily(lat, lon, w$start_ratoon, w$end_ratoon)

    tibble(
      data_id = data_id,
      main_t2m = main_clim$t2m,
      main_ppt = main_clim$ppt,
      main_srad = main_clim$srad,
      ratoon_t2m = ratoon_clim$t2m,
      ratoon_ppt = ratoon_clim$ppt,
      ratoon_srad = ratoon_clim$srad
    )
  })

  saveRDS(out, climate_cache)
  out
}

climate <- get_climate_data(raw)

#-----------------------------
# 4. 合并 + 分类变量
#-----------------------------
dat <- raw |>
  left_join(soil, by = c("lat", "lon")) |>
  left_join(climate, by = "data_id") |>
  mutate(
    clay_class = cut(
      soil_clay_gkg,
      breaks = c(-Inf, 200, 350, Inf),
      labels = c("sandy(<20%)", "loam(20-35%)", "clayey(>35%)")
    ),
    ph_class = cut(
      soil_ph,
      breaks = c(-Inf, 5.5, 7.5, Inf),
      labels = c("acid", "neutral", "alkaline")
    ),
    som_class = cut(
      soil_som_gkg,
      breaks = c(-Inf, 20, 40, Inf),
      labels = c("low", "medium", "high")
    ),
    n_total_class = cut(
      n_total,
      breaks = c(-Inf, 150, 250, Inf),
      labels = c("low_N", "medium_N", "high_N")
    ),
    main_temp_class = cut(
      main_t2m,
      breaks = c(-Inf, 22, 28, Inf),
      labels = c("cool", "optimum", "hot")
    ),
    ratoon_temp_class = cut(
      ratoon_t2m,
      breaks = c(-Inf, 20, 26, Inf),
      labels = c("cool", "optimum", "hot")
    )
  )

write_csv(dat, "results/merged_dataset.csv")

#-----------------------------
# 5. 探索性分析
#-----------------------------
missing_tbl <- dat |>
  summarise(across(everything(), ~sum(is.na(.)))) |>
  pivot_longer(cols = everything(), names_to = "variable", values_to = "n_missing") |>
  mutate(p_missing = n_missing / nrow(dat)) |>
  arrange(desc(p_missing))

write_csv(missing_tbl, "results/missing_summary.csv")

target_vars <- c("m_yield", "r_yield", "m_panicle", "m_spikelet", "r_panicle", "r_spikelet")
normality <- map_dfr(target_vars, function(v) {
  x <- dat[[v]]
  x <- x[!is.na(x)]
  if (length(x) < 3) {
    return(tibble(variable = v, n = length(x), shapiro_p = NA_real_))
  }
  p <- tryCatch(shapiro.test(sample(x, min(length(x), 5000)))$p.value, error = function(e) NA_real_)
  tibble(variable = v, n = length(x), shapiro_p = p)
})

write_csv(normality, "results/normality_shapiro.csv")

#-----------------------------
# 6. 建模：解释再生季产量
#-----------------------------
model_df <- dat |>
  select(
    r_yield, m_yield, m_panicle, m_spikelet, m_grain_wt, m_filled_rate,
    r_panicle, r_spikelet, r_grain_wt, r_filled_rate,
    genotype, variety, cut_ht, n_main_crop, n_ratoon, n_total,
    soil_clay_gkg, soil_n_gkg, soil_som_gkg, soil_ph,
    main_t2m, main_ppt, ratoon_t2m, ratoon_ppt,
    country, year
  ) |>
  mutate(year = as.factor(year)) |>
  drop_na(r_yield)

# 6.1 LMM
lmm_formula <- r_yield ~ m_yield + m_panicle + m_spikelet + m_grain_wt + m_filled_rate +
  cut_ht + n_main_crop + n_ratoon +
  soil_clay_gkg + soil_n_gkg + soil_som_gkg + soil_ph +
  main_t2m + main_ppt + ratoon_t2m + ratoon_ppt +
  (1 | variety) + (1 | country) + (1 | year)

lmm_fit <- lmer(lmm_formula, data = model_df, REML = FALSE)
write_csv(broom.mixed::tidy(lmm_fit, effects = "fixed"), "results/lmm_fixed_effects.csv")
write_csv(as.data.frame(performance::check_collinearity(lmm_fit)), "results/lmm_vif.csv")

# 6.2 鲁棒回归（对异常值敏感性分析）
robust_fit <- MASS::rlm(
  r_yield ~ m_yield + m_panicle + m_spikelet + m_grain_wt + m_filled_rate +
    cut_ht + n_main_crop + n_ratoon +
    soil_clay_gkg + soil_n_gkg + soil_som_gkg + soil_ph +
    main_t2m + main_ppt + ratoon_t2m + ratoon_ppt,
  data = model_df
)
robust_coef <- summary(robust_fit)$coefficients |>
  as.data.frame() |>
  rownames_to_column("term")
write_csv(robust_coef, "results/robust_rlm_coefficients.csv")

# 6.3 LASSO筛选
x <- model.matrix(
  r_yield ~ m_yield + m_panicle + m_spikelet + m_grain_wt + m_filled_rate +
    r_panicle + r_spikelet + r_grain_wt + r_filled_rate +
    genotype + cut_ht + n_main_crop + n_ratoon + n_total +
    soil_clay_gkg + soil_n_gkg + soil_som_gkg + soil_ph +
    main_t2m + main_ppt + ratoon_t2m + ratoon_ppt,
  data = model_df
)[, -1]
y <- model_df$r_yield

cvfit <- cv.glmnet(x, y, alpha = 1, nfolds = 10, standardize = TRUE)
lasso_coef <- coef(cvfit, s = "lambda.1se") |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column("term") |>
  rename(coef = `1`) |>
  filter(coef != 0)
write_csv(lasso_coef, "results/lasso_selected_coefficients.csv")

# 6.4 随机森林
rf_df <- model_df |>
  select(-country, -year, -variety) |>
  drop_na()

rf_fit <- ranger(
  r_yield ~ ., data = rf_df,
  importance = "permutation",
  num.trees = 1000,
  seed = 2026
)
rf_imp <- enframe(rf_fit$variable.importance, name = "variable", value = "importance") |>
  arrange(desc(importance))
write_csv(rf_imp, "results/rf_variable_importance.csv")

#-----------------------------
# 7. 图形输出
#-----------------------------
p1 <- ggplot(dat, aes(x = m_yield, y = r_yield)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "头季产量", y = "再生季产量", title = "头季与再生季产量关系") +
  theme_bw(base_size = 12)

ggsave("results/figures/m_vs_r_yield_scatter.png", p1, width = 7, height = 5, dpi = 300)

p2 <- missing_tbl |>
  mutate(variable = fct_reorder(variable, p_missing)) |>
  ggplot(aes(x = variable, y = p_missing)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "缺失比例", title = "变量缺失情况") +
  theme_bw(base_size = 11)

ggsave("results/figures/missingness_bar.png", p2, width = 8, height = 8, dpi = 300)

message("分析完成。结果已输出到 results/ 目录。")
