## To extract projection temperature & uncertainty from climate models in cmip6 data ##
## Contain 2018-2100 climate model output and Uncertainty ##
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
parallel <- T # parallel switch
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
#Create SpatVector from POIs
crdref <- "+proj=longlat +datum=WGS84"                #坐标参照系统为WGS84
beijing_poi_SpV <- vect(beijing_poi, geom=c("lon", "lat"), crs = crdref)

#fitter "r1i1p1f1" ensample model
tProDir <- str_subset(dir(path_monthly_tPro),"r1i1p1f1") #18 models lefts
##### Parallel setting #####
if(parallel){
  library(foreach)
  library(parallel)
  library(doSNOW)
  myCluster <- makeCluster(4)
  registerDoSNOW(myCluster)
}
##### Main process #####
#load climate_mp data(~1.5mins)


recycle_list <- c(2018:2100)
# if(parallel){
#   iterations <- length(recycle_list)
#   pb <- txtProgressBar(max = iterations, style = 3) # create progress bar
#   progress <- function(n) setTxtProgressBar(pb, n) # update progress bar
#   opts <- list(progress = progress)
# }
pb <- progress_bar$new(total = length(recycle_list)/4, format = "[:bar] :percent")

tic()
## extract delta T by beijing POI in year n
mclapply(recycle_list,function(n){
tProS <- lapply(1:18, function(i){
  nc_monthly_tPro <- rast(paste0(path_monthly_tPro, tProDir[i]))
  tProS <- nc_monthly_tPro
  # Differences
  tProS[[c(((n-2017)*12+1):((n-2017)*12+12))]] = tProS[[c(((n-2017)*12+1):((n-2017)*12+12))]] - tProS[[c(1:12)]]
  tProS <- tProS[[((n-2017)*12+1):((n-2017)*12+12)]]
  for (j in 1:12) names(tProS[[j]]) = j
  tProS<- tProS %>%
    terra::extract(beijing_poi_SpV, bind = T) %>%
    terra::as.data.frame() %>%
    mutate(model = dir(path_monthly_tPro)[i] %>% str_split("_", simplify = T) %>% .[,3]) %>% 
    pivot_longer(cols = c(2:13), names_to = "month", values_to = "deltaT") %>%
    select(id, month, deltaT, model)
  return(tProS)
}) %>% data.table::rbindlist()

#calculate climate model uncerrtainty
### 
tProS <- lapply(1:12, function(i){
  deltaT <- tProS[month == i]
  deltaT[month == i,c(paste0("deltaT",n,".sd"),paste0("deltaT",n,".mean")) := bootstrapCI(deltaT),by = .(id)]
}) %>% data.table::rbindlist()
## clean up deltaT
tProS[, deltaT := NULL]
fwrite(tProS,file = paste0(odir,"tProCI_",n,".csv.gz"), compress = "gzip")
# progress(n)
gc()
pb$tick()
},mc.cores = 4)

toc()
##### Main process 2 #####
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
  
  #筛选2050, 2100
  tProM <- tProM[[c(((2050-2017)*12+1):((2050-2017)*12+12), ((2100-2017)*12+1):((2100-2017)*12+12))]]
  
  #Rename
  for (j in 1:12) {
    names(tProM[[j]]) = paste0("mid_",j)
    names(tProM[[j+12]]) = paste0("long_",j)
  }
  tProM<- tProM %>%
    terra::extract(beijing_poi_SpV, bind = T) %>%
    as.data.frame() %>%
    mutate(model = dir(path_monthly_tPro)[i])
  return(tProM)
}, mc.cores = 4) %>% data.table::rbindlist()
toc()
# extract model names
tProM[,model := model %>% str_split("_", simplify = T) %>% .[,3]]

mid <- tProM %>%
  pivot_longer(cols = c(2:13), names_to = "month", values_to = "deltaTMid") %>%
  separate(col = c("month"), sep = "_", into = c("trash", "month"), convert = T) %>%
  select(id, month, deltaTMid, model)
long <- tProM %>%
  pivot_longer(cols = c(14:25), names_to = "month", values_to = "deltaTLong") %>%
  separate(col = c("month"), sep = "_", into = c("trash", "month"), convert = T) %>%
  select(id, month, deltaTLong, model)

tProM <- left_join(mid,long,by = c("id", "month", "model"))

rm(list = c("mid", "long"))
tProM <- tProM %>% as.data.table()
# Bootstrapping
tProM <- mclapply(1:12, function(i){
  deltaT <- tProM[month == i]
  deltaT[month == i,c("deltaTMid.sd","deltaTMid.mean") := bootstrapCI(deltaTMid),by = .(id)]
  deltaT[month == i,c("deltaTLong.sd","deltaTLong.mean") := bootstrapCI(deltaTLong),by = .(id)]
  return(deltaT)
}, mc.cores = 6) %>% data.table::rbindlist()
tProM <- tProM[,.(id, month, deltaTMid.sd, deltaTMid.mean, deltaTLong.sd, deltaTLong.mean)]

# Join 2018-2100 timeline and 21st centrary mid|long
# tProM <- fread("~/ClimBehav/results/data/processed/tProMax.csv.gz",nThread = 14)
#clean up tProM
tPro <- tProM %>% distinct(id, month, .keep_all = T)

for (n in 2018:2100) {
  deltaT <- fread(paste0(odir,"tProCI_",n,".csv.gz"),nThread = 14) %>% 
    distinct(id, month, .keep_all = T)
  deltaT <- deltaT[,model := NULL]
  tPro <- left_join(tPro, deltaT,by = c("id","month"))
}
## fwrite tPro in year n or mid|long
fwrite(tPro, file = "~/ClimBehav_bj/results/data/processed/tProCI.csv.gz", compress = "gzip", nThread = 14)
##### Stop parallel #####
# if(parallel) stopCluster(myCluster)