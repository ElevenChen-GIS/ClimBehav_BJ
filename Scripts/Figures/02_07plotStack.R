#Stacking diagrame
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
library(cowplot) #拼图 for insert_xaxis_grob
library(scales) # for scales axis
load("~/ClimBehav_bj/results/data/raw/01_AOI/11-北京市-aoi.Rdata") #load AOI

threshold <- function(dt, a){
  visit.stat <- dt %>% group_by(id) %>% dplyr::summarise(daily_visit_num = sum(visits / 365, na.rm = T))
  base_select <- filter(visit.stat, daily_visit_num >= a)
  dt <- dt %>%
    dplyr::filter(id %in% base_select$id)
  return(dt)
}
squash_axis <- function(from, to, factor) { 
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}

#读取climate_mp数据
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
gc()

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

#######4.Work#########
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

#手机位置请求百分比堆积图
visit.stat <- df %>% 
  group_by(class, date) %>% dplyr::summarise(visit_num = sum(visits, na.rm = T)) %>% 
  setDT()
visit.stat[, visit_all := sum(visit_num, na.rm = T), by = date]
visit.stat[, visit_perc := visit_num/visit_all]

#transfer Transport to 1st position

visit.stat$class <- visit.stat$class %>% 
  factor(levels = c("Transport", "Residential", "Park", "Work"))
pStack <- ggplot(visit.stat) + 
  geom_area(aes(x = date, y = visit_perc, fill = class),alpha = 0.7) +
  scale_fill_manual(values = c("Residential" = "#457b9d","Park" = "#84a729","Transport" = "#C75c64","Work" = "#f0b57d"))+
  labs(y = "Persentage of visitation",
       x = "Date",
       fill = "Class\nof POIs") +
  theme_classic()+
  scale_x_date(date_breaks = "1 month",date_labels =c("","J","F","M","A","M","J","J","A","S","O","N","D","J"))+
  guides(fill = guide_legend(theme = theme(
    text = element_text(size = 15)
  ))) +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12))+
  coord_trans(y = squash_axis(0.3, 1, 30))
save(pStack, file = "~/ClimBehav_bj/results/data/processed/pStack.RData")
