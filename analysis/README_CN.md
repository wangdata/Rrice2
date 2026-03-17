# 再生稻 Meta 数据分析说明（R）

## 你将得到什么
运行 `analysis/ratoon_analysis.R` 后，会在 `results/` 生成：

- `merged_dataset.csv`：合并原始数据 + 土壤 + 气象 + 分类变量。
- `missing_summary.csv`：每列缺失值统计。
- `normality_shapiro.csv`：目标变量正态性检验（Shapiro-Wilk）。
- `lmm_fixed_effects.csv`：线性混合模型固定效应结果。
- `lmm_vif.csv`：共线性检查结果。
- `robust_rlm_coefficients.csv`：鲁棒回归系数。
- `lasso_selected_coefficients.csv`：LASSO筛选变量。
- `rf_variable_importance.csv`：随机森林变量重要性。
- `figures/`：散点图与缺失图。

## 运行步骤
1. 安装 R（建议 >= 4.3）。
2. 安装依赖包：

```r
install.packages(c(
  "tidyverse", "readxl", "janitor", "httr2", "jsonlite", "lubridate",
  "nasapower", "lme4", "lmerTest", "MASS", "glmnet", "ranger",
  "performance", "broom.mixed", "car"
))
```

3. 在仓库根目录运行：

```bash
Rscript analysis/ratoon_analysis.R
```

## 关于“HWSD”说明
脚本默认通过 SoilGrids REST API 按经纬度抓取土壤属性（黏粒、氮、有机碳、pH），这是公开可用且自动化程度高的方式。

如果你坚持使用 HWSD 最新栅格，请在脚本中将 `query_soilgrids()` 部分替换为：
- 本地下载 HWSD 对应变量栅格（clay、N、SOC/OM、pH）；
- 使用 `terra::extract()` 对每个样点提取数值；
- 按脚本相同字段名合并。

## 头季/再生季气象窗口说明
当前数据缺少精确播栽/收获日期，脚本用如下近似规则：
- 北半球头季起点：4月1日；南半球：10月1日；
- 头季长度采用 `Mgrowth length`；
- 再生季紧接头季结束后，以 `Rgrowth length` 为长度。

建议你后续补充“真实日期列”替换该规则，以提升气象匹配精度。

