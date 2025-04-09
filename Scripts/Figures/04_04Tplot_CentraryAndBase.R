##### Initialization #####
setwd("/project/public/tmp_cy/")
library(tidyverse)
library(terra)       #处理栅格数据
library(parallel)     #提供并行运算
library(sf)
library(dplyr)
library(data.table)
library(tictoc)
library(progress)
library(plyr) #for additional manipulation functions and revulue function
bootstrapCI <- function(x){
  mean_samples <- replicate(1000, expr = {
    y = sample(x, 18,replace = T)
    mean(y)
  })
  a <- sd(mean_samples)/sqrt(length(mean_samples))
  b <- mean(mean_samples)
  return(list(a,b))
}
crs_84 <- st_crs("EPSG:4326")  ## WGS 84 大地坐标
cores <- 4 #paralell in 4 cores
parallel <- T
odir <- "tProYears/"  # dir of results
path_monthly_tPro <- "/project/public/tmp_cy/tProModels/"

#导入Beijing POI
load("~/ClimBehav_bj/results/data/raw/01_AOI/11-北京市-aoi.Rdata")

# Classification and filter
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
aoi <- inner_join(df_classify, aoi, by = "id", relationship = "many-to-many")

#sample 10 pois per class, then a 40 pois table
beijing_poi <- rbind(
  aoi[id %in% sample(aoi[class == 1]$id,10)],
  aoi[id %in% sample(aoi[class == 2]$id,10)],
  aoi[id %in% sample(aoi[class == 3]$id,10)],
  aoi[id %in% sample(aoi[class == 4]$id,10)]
)

beijing_poi <- beijing_poi %>%
  as_tibble() %>%
  separate(col = "location", into = c("lon", "lat"), sep = ",", convert = T) %>%
  select(id, lon, lat) %>%
  distinct(id, .keep_all = T)                          #剔除重复的id,提出后为40条观测,10 per class
rm("aoi")
#load temperature data in 2017
# tic()
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
df <- df %>% dplyr::filter(hour >= 7 & hour <= 22)
#extract temperature by 40 pois
beijingT_poi <- left_join(beijing_poi,df,by = "id")
rm(df)
gc()

#Create SpatVector from POIs
crdref <- "+proj=longlat +datum=WGS84"                #坐标参照系统为WGS84
beijing_poi_SpV <- vect(beijing_poi, geom=c("lon", "lat"), crs = crdref)

#fitter "r1i1p1f1" ensample model
tProDir <- str_subset(dir(path_monthly_tPro),"r1i1p1f1") #18 models lefts

##### Main process #####
recycle_list <- c(1:18)

tic()
## 21st century monthly pre(~ 20s)
tProM <- mclapply(1:18, function(i){
  nc_monthly_tPro <- rast(paste0(path_monthly_tPro, tProDir[i]))
  ##Average by month
  tProM <- NULL
  tProM <- nc_monthly_tPro
  for (j in 1:12) {
    tProM[[(2050-2017)*12+j]] = mean(nc_monthly_tPro[[(2046:2055-2017)*12+j]])
    tProM[[(2100-2017)*12+j]] = mean(nc_monthly_tPro[[(2090:2099-2017)*12+j]])
  }
  
  ## Differences
  for (n in c(2050, 2100)){
    for (j in 1:12) {
      tProM[[(n-2017)*12+j]] = tProM[[(n-2017)*12+j]] - tProM[[j]]
    }
  }
  
  #筛选2017, 2050, 2100
  tProM <- tProM[[c(1:12, ((2050-2017)*12+1):((2050-2017)*12+12), ((2100-2017)*12+1):((2100-2017)*12+12))]]
  
  #Rename
  for (j in 1:12) {
    names(tProM[[j]]) = paste0("base_",j)
    names(tProM[[j+12]]) = paste0("mid_",j)
    names(tProM[[j+24]]) = paste0("long_",j)
  }
  tProM<- tProM %>%
    terra::extract(beijing_poi_SpV, bind = T) %>%
    as.data.frame() %>%
    mutate(model = dir(path_monthly_tPro)[i])
  return(tProM)
}, mc.cores = cores) %>% data.table::rbindlist()
toc()
# extract model names
tProM[,model := model %>% str_split("_", simplify = T) %>% .[,3]]


base <- tProM %>%
  pivot_longer(cols = c(2:13), names_to = "month", values_to = "Tbase") %>%
  separate(col = c("month"), sep = "_", into = c("trash", "month"), convert = T) %>%
  select(id, month, Tbase, model)
mid <- tProM %>%
  pivot_longer(cols = c(14:25), names_to = "month", values_to = "deltaTMid") %>%
  separate(col = c("month"), sep = "_", into = c("trash", "month"), convert = T) %>%
  select(id, month, deltaTMid, model)
long <- tProM %>%
  pivot_longer(cols = c(26:37), names_to = "month", values_to = "deltaTLong") %>%
  separate(col = c("month"), sep = "_", into = c("trash", "month"), convert = T) %>%
  select(id, month, deltaTLong, model)

tProM <- left_join(base, mid,by = c("id", "month", "model")) %>% left_join(long,by = c("id", "month", "model"))

rm(list = c("mid", "long","base"))
tProM <- inner_join(beijingT_poi,tProM,by = c("id","month"), relationship = "many-to-many")
tProM <- tProM %>% as.data.table()
# Bootstrapping
tProM[,Tbase := t]
tProM[,TMid := Tbase + deltaTMid]
tProM[,TLong := Tbase + deltaTLong]
tProM[,':='(
  Tbase.mean = min(Tbase),
  TMid.mean = min(TMid),
  TLong.mean = min(TLong)
),by = .(month)]
tProM <- tProM %>% distinct(month,.keep_all = T)
#remove year and t colums
tProM$year <- NULL
tProM$t <- NULL
tProM <- tProM %>% 
  pivot_longer(cols = c("Tbase.mean","TMid.mean","TLong.mean"), names_to = "year", values_to = "t") %>%
  select(month, t, year)

tProM$year <- tProM$year %>% 
  factor() %>%
  revalue(c("Tbase.mean" = "2017", "TMid.mean" = "2046-2055", "TLong.mean" = "2090-2099")) %>%
  factor(levels = c("2017","2046-2055", "2090-2099"))
# tPro[,date := make_date(year = 2017, month, day)]
# tPro[,doy := as.integer(date-make_date(year = 2017,month = 1,day =1))+1]
p <- ggplot(tProM,aes(month,t))+
  geom_line(aes(color = year))+
  scale_colour_manual(values = c("2017"="#457b9d","2046-2055"="#84a729", "2090-2099"="#C75c64"))+
  scale_x_continuous(breaks = c(1:12)) +
  ylab("Temperature(℃)")+
  theme_bw()
ggsave("~/ClimBehav_bj/results/plots/Tpro.pdf",plot = p, device = cairo_pdf, width =6, height =6)
ggsave("~/ClimBehav_bj/results/plots/Tpro.eps",plot = p, device = eps, width = 6, height = 6)
