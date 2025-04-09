## for Descriptive analysis ##
## created on 2025-03-27    ##
## update on 2025-03-27     ##

---
  
#Set up packages and environment
library(tidyverse)
library(data.table)
load("~/ClimBehav_bj/results/data/raw/01_AOI/11-北京市-aoi.Rdata")
threshold <- function(dt, a){
  visit.stat <- dt %>% group_by(id) %>% dplyr::summarise(daily_visit_num = sum(visits / 365, na.rm = T))
  base_select <- filter(visit.stat, daily_visit_num >= a)
  dt <- dt %>%
    dplyr::filter(id %in% base_select$id)
  return(dt)
}
#--------------------------------------描述性统计
mystats <- function(x,na.omit=T){
  if(na.omit)
    x <- x[!is.na(x)]
  n <- length(x)
  mean <- mean(x)
  sd <- sd(x)
  min <- min(x)
  qu_1st <- stats::quantile(x, probs = c(0.25))
  qu_3st <- stats::quantile(x, probs = c(0.75))
  max <- max(x)
  return(c(n = n, mean = mean, sd = sd, min = min, qu_1st, qu_3st,max =  max))
}

#Data preperation 
## Climate & MP
#读取climate_mp数据(~1.5mins)
df <- fread(file = "~/ClimBehav_bj/results/data/processed/climate_mp.csv.gz", nThread = 14)

#df数据集还原
dig = 3
cutt_levels <- c("(-10,-5]", "(-5,0]", "(0,5]", "(5,10]", "(10,15]", "(15,20]", "(20,25]", "(25,30]", "(30,35]", "(35,40]")

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
  cutt = factor(cut(t, breaks = c(seq(from = -10, to = 40, by = 5))), levels = cutt_levels))]
# toc()

##Classification 
# 这一部分主要为了探究如何从mp数据中找到能代表全域的AOI，并且将现有类别进行clear
# 1.  classify AOI using behavior system;
# 
# 2.  filter out large AOI with visitors greater than 50 people per day;
# 
# 3.  filter out the daytime(7 am ~ 23pm).

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
gc()

#threshold by class
df <- df %>% threshold(50) %>%
  dplyr::filter(hour >= 7 & hour <= 22)
gc()

#statistic df
stat_df <- sapply(df[,.(visits, t, aqi, tp, wind, rhum, sp, tcc, blh)], mystats)

#save as csv
write.csv(stat_df, file = "~/ClimBehav_bj/results/plots/statDf.csv")
