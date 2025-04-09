###这里需要先运行02_03classRegression中的基础（大循环前的）###
# library(progress)

#随机抽样参数
set.seed(123)  # 设置随机种子以便复现
sampling_size <- 1000

residual_plot <- NULL

## 第一个位置：新建一个其实进度条
pb <- txtProgressBar(style=3)
recycle_list <- c(1:4)

star_time <- Sys.time() ## 记录程序开始时间

for(i in 1:4){
    x <- workLiveLabel[i]
    dt <- df[class == x]
    dt <- within(dt, cutt <- relevel(cutt, ref = "(20,25]"))
    
    reg[[i]] <- felm(log_visit ~
                       cutt + tp + wind + rhum + sp + tcc + holidays +aqi +blh        #   covariates |
                     |id + hour + month + dow
                     |0|type,
                     data = dt, na.action="na.omit")
  fitted_values <- fitted(reg[[i]])  # 提取拟合值
  residuals <- resid(reg[[i]])  # 提取残差
  
  # 随机抽样（例如，抽取1000个观测值）
  sample_indices <- sample(1:nrow(residuals), sampling_size)
  fitted_values_sample <- fitted_values[sample_indices, ]
  residuals_sample <- residuals[sample_indices, ]
  residual_plot[[i]] <- ggplot(data = data.frame(fitted_values_sample, residuals_sample), aes(x = fitted_values_sample, y = residuals_sample)) +
    geom_point() +  # 添加散点
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # 添加参考线
    labs(title = paste0(x,":Residuals vs Fitted Values"),
         x = "Fitted Values",
         y = "Residuals") +
    theme_minimal()
  ## 第二个位置：实时反映进度
  setTxtProgressBar(pb, i/length(recycle_list))
}

end_time <- Sys.time()  ## 记录程序结束时间

## 第三个位置关闭进度条
close(pb)

run_time <- end_time - star_time  ## 计算程序运行时间

library(stargazer)
stargazer(reg[[1]],reg[[2]],reg[[3]],reg[[4]],
          # coef_holiday_list[[2]], coef_holiday_list[[3]], coef_holiday_list[[4]],
          dep.var.labels.include = F,
          # omit.stat = c("LL","ser","f","rsq"),
          keep = c(1:9),
          star.cutoffs = c(0.1, 0.05, 0.01),
          type = "text"
)
