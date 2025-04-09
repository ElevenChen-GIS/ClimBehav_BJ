library(tidyverse)
# library(terra)       #处理栅格数据
library(magrittr)    #提供管道操作
library(lubridate)    #提供时间有关函数
library(data.table)
library(fixest)   #提供feols函数
library(lfe)  #提供felm函数
# library(ggExtra) # for marginal histograms
library(plyr) #for additional manipulation functions
# library(gridExtra) # for multiplots
# library(tictoc) # measure time
# library(parallel)     #提供并行运算
library(zoo) #提供滑动平均函数
# library(grid) #自定义常用主题
# library(cowplot) #拼图
library(stargazer)

source("~/ClimBehav_bj/results/scripts/binned_plot_function.R")
load("~/ClimBehav_bj/results/data/raw/01_AOI/11-北京市-aoi.Rdata")
threshold <- function(dt, a){
  visit.stat <- dt %>% group_by(id) %>% dplyr::summarise(daily_visit_num = sum(visits / 365, na.rm = T))
  base_select <- filter(visit.stat, daily_visit_num >= a)
  dt <- dt %>%
    dplyr::filter(id %in% base_select$id)
  return(dt)
}
bootstrapCI <- function(n,a,b){
  y = NULL
  delta <- replicate(1000,expr = {
    rnorm(n, a, b)})
  for (i in 1:nrow(delta)) {
    y <- rbind(y, as.numeric(quantile(delta[i,],c(.025,.975))))
  }
  return(list(y[,1],y[,2]))
}
season<- function(df){
  df$season<-4
  df$season[df$month >= 9 & df$month <= 11]<-3
  df$season[df$month >= 6 & df$month <= 8]<-2
  df$season[df$month >= 3 & df$month <= 5]<-1
  return(df)
}
# source("~/ClimBehav_bj/results/scripts/change_plot.R")

#读取climate_mp数据(~1.5mins)
df <- fread(file = "~/ClimBehav_bj/results/data/processed/climate_mp.csv.gz", nThread = 14)

#df数据集还原
dig = 3
cutt_levels <- c("(-10,-5]", "(-5,0]", "(0,5]", "(5,10]", "(10,15]", "(15,20]", "(20,25]", "(25,30]", "(30,35]", "(35,40]")
cutt_levels2 <- c("(-9,-6]", "(-6,-3]", "(-3,0]", "(0,3]", "(3,6]", "(6,9]", "(9,12]", "(12,15]", "(15,18]","(18,21]", "(21,24]", "(24,27]", "(27,30]", "(30,33]", "(33,36]", "(36,39]") 

df[,':='(
  date = make_date(year, month, day),
  t = t / 10 ^ dig,
  sp = sp / 10^dig,   #将单位换算成hPa
  tp = tp * 10^3,    #将单位换算成mm
  blh = blh / 10^dig, 
  tcc = tcc / 10^3,     #将单位换算成%   
  wind = wind / 10^dig,
  rhum = rhum / 10^dig,
  visits = as.double(visits)
)][, ':='(
  log_visit = log(visits + 1),
  cutt = factor(cut(t, breaks = c(seq(from = -10, to = 40, by = 5))), levels = cutt_levels),
  cutt2 = factor(cut(t, breaks = c(seq(from = -9, to = 39, by = 3))), levels = cutt_levels2))]


###### 1.classification ######
##############################

df_classify<- aoi %>% as.data.table()
#######1.居家行为#########
df_classify[str_detect(type, "^120302") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 1, subclass = 1)] #住宅小区

#######2.公园访问#########
df_classify[str_detect(type, "^1101(00|01|05|06)") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 2, subclass = 2)] #公园广场和内部设施
df_classify[str_detect(type, "^1101(02|03|04)") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 2, subclass = 3)] #主题公园
df_classify[str_detect(type, "^1102(00|01|02|03|04|05|06|07|09)") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 2, subclass = 4)] #风景名胜

####### 3.Transport #########
df_classify[str_detect(type, "^(150104|150200|1504)") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 3, subclass = 5)] #机场、火车站、长途车站
df_classify[str_detect(type, "^1503(01|02|03)|151200") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 3, subclass = 6)] #渡口、港口

#######4.Workplace#########
df_classify[str_detect(type, "^1704") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 4, subclass = 7)] #农林牧渔基地
df_classify[str_detect(type, "^1202|^17(01|02)") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 4, subclass = 8)] #公司企业、写字楼
df_classify[str_detect(type, "^1703|^1201") & str_detect(name, "(建设中)|(装修中)") == 0, ':='(class = 4, subclass = 9)] #工厂，产业园
df_classify[str_detect(name, "(建设中)|(装修中)"), ':='(class = 4, subclass = 10)] #建筑工地,需要放在所有分类的最后
df_classify <- df_classify[is.na(class) == 0, .(id, class, subclass)]
df <- inner_join(df_classify, df, by = "id", relationship = "many-to-many")
rm(df_classify)
#threshold by class
df <- df %>% threshold(50) %>%
  dplyr::filter(hour >= 7 & hour <= 22)
gc()


reg <- NULL
reg_date <- NULL
reg_season <- NULL
reg_3bin <- NULL

t_coef <- NULL
t_coef_date <- NULL
t_coef_season <- NULL
t_coef_3bin <- NULL

workLiveLabel <- c("Residential", "Park", "Transport", "Work")
df$class <- df$class %>% 
  factor() %>% 
  revalue(c("1" = "Residential", "2" = "Park", "3" = "Transport", "4" = "Work")) %>% 
  factor(levels = c("Residential", "Park", "Transport", "Work"))
df <- season(df)
## 第一个位置：新建一个其实进度条
pb <- txtProgressBar(style=3)
recycle_list <- c(1:4)

star_time <- Sys.time() ## 记录程序开始时间

for (i in 1:4) {
  x <- workLiveLabel[i]
  dt <- df[class == x]
  dt <- within(dt, cutt <- relevel(cutt, ref = "(20,25]"))
  
  #base
  reg[[i]] <- felm(log_visit ~
                     cutt + tp + wind + rhum + sp + tcc + holidays +aqi +blh
                   |id + hour + month + dow
                   |0|id,
                   data = dt, na.action="na.omit")
  # summary(reg[[i]])
  t_coef[[i]] <- felmCoeff(reg[[i]],var = "cutt", omit = c(20, 25)) %>%
    select(coef, se, star, xmid, range, ci.l, ci.h)
  deltaReg <- as.data.table(t_coef[[i]])
  deltaReg[,c("ci.l","ci.h") := bootstrapCI(nrow(deltaReg),coef, se)]
  t_coef[[i]] <- deltaReg
  rm(deltaReg)
  gc()
  
  #date FE instead of DOW FE
  reg_date[[i]] <- felm(log_visit ~
                     cutt + tp + wind + rhum + sp + tcc + holidays +aqi +blh
                   |id + hour + month + date
                   |0|id,
                   data = dt, na.action="na.omit")
  # summary(reg[[i]])
  t_coef_date[[i]] <- felmCoeff(reg_date[[i]],var = "cutt", omit = c(20, 25)) %>%
    select(coef, se, star, xmid, range, ci.l, ci.h)
  deltaReg <- as.data.table(t_coef_date[[i]])
  deltaReg[,c("ci.l","ci.h") := bootstrapCI(nrow(deltaReg),coef, se)]
  t_coef_date[[i]] <- deltaReg
  rm(deltaReg)
  gc()
  
  #season instead of month FE
  reg_season[[i]] <- felm(log_visit ~
                          cutt + tp + wind + rhum + sp + tcc + holidays +aqi +blh
                        |id + hour + season + dow
                        |0|id,
                        data = dt, na.action="na.omit")
  # summary(reg[[i]])
  t_coef_season[[i]] <- felmCoeff(reg_season[[i]],var = "cutt", omit = c(20, 25)) %>%
    select(coef, se, star, xmid, range, ci.l, ci.h)
  deltaReg <- as.data.table(t_coef_season[[i]])
  deltaReg[,c("ci.l","ci.h") := bootstrapCI(nrow(deltaReg),coef, se)]
  t_coef_season[[i]] <- deltaReg
  rm(deltaReg)
  gc()
  
  # 3 bins instead of 5 bins
  dt <- within(dt, cutt2 <- relevel(cutt2, ref = "(21,24]"))
  reg_3bin[[i]] <- felm(log_visit ~
                            cutt2 + tp + wind + rhum + sp + tcc + holidays +aqi +blh
                          |id + hour + month + dow
                          |0|id,
                          data = dt, na.action="na.omit")
  # summary(reg[[i]])
  t_coef_3bin[[i]] <- felmCoeff(reg_3bin[[i]],var = "cutt2", omit = c(21, 24)) %>%
    select(coef, se, star, xmid, range, ci.l, ci.h)
  deltaReg <- as.data.table(t_coef_3bin[[i]])
  deltaReg[,c("ci.l","ci.h") := bootstrapCI(nrow(deltaReg),coef, se)]
  t_coef_3bin[[i]] <- deltaReg
  rm(deltaReg)
  gc()
  setTxtProgressBar(pb, i/length(recycle_list))
}

end_time <- Sys.time()  ## 记录程序结束时间
## 第三个位置关闭进度条
close(pb)

run_time <- end_time - star_time  ## 计算程序运行时间

save(t_coef,t_coef_date,t_coef_season,t_coef_3bin, file = "~/ClimBehav_bj/results/data/processed/pClassRobust.RData")