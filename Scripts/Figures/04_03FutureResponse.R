##To project future response change and uncertainty from 2017 to 2018-2100 or mid|long 21st centrary##
##Use data from 2018-2100 climate model output and Uncertainty                                      ##
##### Initialization #####
setwd("/project/public/tmp_cy/")
library(tidyverse)
library(tictoc)
library(parallel)     #提供并行运算
library(sf)
library(sfheaders) #for sf_point function
library(dplyr)
library(data.table)
library(splines)      #提供spline function
library(lubridate)    #提供时间有关函数
library(progress) #provide progress bar function
library(RANN)

#filter daily visitation of poi
threshold <- function(dt, a){
  visit.stat <- dt %>% group_by(id) %>% dplyr::summarise(daily_visit_num = sum(visits / 365, na.rm = T))
  base_select <- filter(visit.stat, daily_visit_num >= a)
  dt <- dt %>%
    dplyr::filter(id %in% base_select$id)
  return(dt)
}
crs_84 <- st_crs("EPSG:4326")  ## WGS 84 大地坐标
cores = 9  #call for 9 cores
parallel <- T # parallel switch
load("~/ClimBehav_bj/results/data/raw/01_AOI/11-北京市-aoi.Rdata") # prepare to subclass
odir <- "respProj/"  # dir of results
rep <- 100
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
df <- fread(file = "~/ClimBehav_bj/results/data/processed/climate_mp.csv.gz")
dig = 3

df[,':='(
  date = make_date(year, month, day),
  t = t / 10 ^ dig,
  visits = as.double(visits)
)][, ':='(
  log_visit = log(visits + 1)
)]
df <- df[,.(id, type, year, month, day, hour, t, visits, date, log_visit)]
rm("dig")

#Load regression coef & monthly projection temperature
tProCI <- fread("~/ClimBehav_bj/results/data/processed/tProCI.csv.gz")
tCoef <- fread("~/ClimBehav_bj/results/data/processed/t_coef.csv")

# extend tProCI
## trans aoi to beijing poi
beijing_poi <- aoi %>%
  as_tibble() %>%
  separate(col = "location", into = c("lon", "lat"), sep = ",", convert = T) %>%
  select(id, lon, lat) %>%
  distinct(id, .keep_all = T) 
## create beijing poi sf
beijing_point<- sf_point(beijing_poi, x = "lon", y = "lat", keep = T)
st_crs(beijing_point) <- crs_84

## create tProCI sample poi sf
tProCI_poi <- tProCI %>% 
  distinct(id) %>% 
  left_join(beijing_poi,by = "id")

tProCI_point<- sf_point(tProCI_poi, x = "lon", y = "lat", keep = T)
st_crs(tProCI_point) <- crs_84

# ggplot(tProCI_point) +
#   geom_sf()

#碰撞检测
# any(st_intersects(tProCI_point, beijing_point, sparse = FALSE))  #TRUE

#最邻近点使用RANN包
p_beijing=st_coordinates(beijing_point) #获取点坐标。由于这个包只返回欧式距离，不是专门针对地理空间开发的，所以在计算前还是先把点坐标转换成投影坐标系为宜。
p_tProCI=st_coordinates(tProCI_point)
d=nn2(p_tProCI, p_beijing, k=1) #k=2即为设定K级最邻近点，因为这个包只支持计算两个点集之间的最邻近点关系，所以只能计算包含自身的2级最邻近点。该函数还有其他参数可以设定搜寻方法、查寻树的类型、搜索半径等。
d_indice=d[[1]]  #返回值包括两个list元素，第一个元素就是最邻近点序号矩阵,第二列就是我们要找的最邻近点序号。
beijing_poi$near <- d_indice[,1]
tProCI_poi$near <- c(1:nrow(tProCI_poi))
tProCI <- left_join(tProCI,tProCI_poi,by = "id")
aoi <- left_join(aoi, beijing_poi,by = "id")
rm(list = c("beijing_poi", "beijing_point","d","d_indice","p_beijing","p_tProCI","tProCI_poi","tProCI_point"))

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
df_classify <- df_classify[is.na(class) == 0, .(id, class, subclass, near)]
df <- inner_join(df_classify, df, by = "id", relationship = "many-to-many")

rm(list = c("df_classify", "aoi"))

df <- df %>% threshold(50) %>%
  dplyr::filter(hour >= 7 & hour <= 22)
gc()
# Fit on natural cubic spline
model <- NULL
model.ci <- NULL
for (i in 1:4) {
    model[[i]] <- lm(coef * 100 ~ splines::ns(xmid, df = 3), data = tCoef[class == i])
}
for (i in 1:4) {
  for (j in 1:rep) {
    tCoef.pro <- tCoef %>% 
      rowwise() %>% 
      mutate(coef = rnorm(1,coef,se)) %>% 
      as.data.table()
    model.ci[[(i-1)*rep+j]] <- lm(coef * 100 ~ splines::ns(xmid, df = 3), data = tCoef.pro[class == i])
    
  }
}
gc() #free used memory

recycle_list <- c(2018:2100,"Mid","Long")
if(parallel){
  iterations <- length(recycle_list)
  pb <- txtProgressBar(max = iterations, style = 3) # create progress bar
  progress <- function(n) setTxtProgressBar(pb, n) # update progress bar
  opts <- list(progress = progress)
}
#(98612 sec,12.4~27.4h)
tic()

foreach(
  n = recycle_list,
  .packages = c("tidyverse","data.table","dplyr"),
  .errorhandling = "pass",
  .options.snow = opts
) %dopar% {
  for (i in names(df)[names(df) %like% "tPro"]) {
    df <- df[, (i) := NULL]
  }
  ## choose climate data in year n and join
  df <- left_join(df, tProCI[,.SD, .SDcols = (names(tProCI) %like% paste0("deltaT", n)) | (names(tProCI) %in% c("near", "month"))],by = c("near", "month"))
  
  setnames(df, old = c(names(df)[names(df) %like% "deltaT"]), new = c("deltaT.sd","deltaT.mean"))
  
  df2 <- df[,.(
    t=mean(t),
    deltaT.mean = mean(deltaT.mean),
    deltaT.sd = mean(deltaT.sd)
  ),by = .(class, month, day, hour)]
  ##calculate project hourly t
  df2[,':='(
    tPro.mean = t + deltaT.mean,
    tPro.sd = deltaT.sd
  )]
  
  ## clean up df
  for (i in names(df)[names(df) %like% "deltaT"]) {
    df <- df[, (i) := NULL]
  }
  
  # predict response gap in year n
  differN.plot <- lapply(1:4, FUN = function(i){
    dt <- df2[class == i]
    forecastNow <- data.table(stats::predict.lm(model[[i]], data.frame(xmid = dt[, (t)])))
    forecastNow.ci <- lapply(1:rep,function(j){data.table(resp = stats::predict.lm(model.ci[[(i-1)*rep+j]], data.frame(xmid = dt[, (t)])))}) %>% bind_cols() %>%  mutate(id = 1:5840) %>% pivot_longer(1:100,names_to = "times",values_to = "resp") %>% as.data.table()
    forecastNow.ci[,lwr := quantile(resp,.025),by = .(id)]
    forecastNow.ci[,upr := quantile(resp,.975),by = .(id)]
    forecastNow.ci <- forecastNow.ci %>% distinct(id,.keep_all = T)
    forecast.mean <- data.table(stats::predict.lm(model[[i]], data.frame(xmid = dt[, (tPro.mean)])))
    forecast.ci <- 
      lapply(1:rep,function(k){
        dt.pro <- dt %>% 
          rowwise() %>% 
          mutate(tPro.mean = rnorm(1,tPro.mean,tPro.sd)) %>% 
          as.data.table()
        resp.pro <- lapply(1:rep,function(j){data.table(resp = stats::predict.lm(model.ci[[(i-1)*rep+j]], data.frame(xmid = dt.pro[, (tPro.mean)])))}) %>% bind_cols() %>%  mutate(id = 1:5840) %>% pivot_longer(1:100,names_to = "times",values_to = "resp") %>% as.data.table()
        return(resp.pro)
        }) %>% data.table::rbindlist()
    forecast.ci[,lwr := quantile(resp,.025),by = .(id)]
    forecast.ci[,upr := quantile(resp,.975),by = .(id)]
    forecast.ci <- forecast.ci %>% distinct(id,.keep_all = T)
    dt$forecastNow.mean <- forecastNow[,1]
    dt$forecastNow.upr <- forecastNow.ci[,upr]
    dt$forecastNow.lwr <- forecastNow.ci[,lwr]
    dt$forecast.mean <- forecast.mean[,1]
    dt$forecast.upr <- forecast.ci[,upr]
    dt$forecast.lwr <- forecast.ci[,lwr]
    
    #calculate gap
    dt[,differ := forecast.mean-forecastNow.mean]
    dt[,differ.upr := forecast.upr-forecastNow.lwr]
    dt[,differ.lwr := forecast.lwr-forecastNow.upr]
    
    ## get mean for each month
    dt <- dt[, .(
      differ = mean(differ, na.rm = T),
      differ.upr = mean(differ.upr, na.rm = T),
      differ.lwr = mean(differ.lwr, na.rm = T)
    ), by = .(month)][,class := i]
    return(dt)
  }) %>% data.table::rbindlist()
  differN.plot[,year := n]
  fwrite(differN.plot,file = paste0("/project/public/tmp_cy/respProj/","respPro_",n,".csv.gz"), compress = "gzip")
  rm(differN.plot)
  gc()
  return(1)
}

toc()

##### Stop parallel #####
if(parallel) stopCluster(myCluster)
