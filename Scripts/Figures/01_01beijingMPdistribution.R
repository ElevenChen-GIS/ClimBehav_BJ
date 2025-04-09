###summarise checkins_raster,beijing mp data tif###
library(ggplot2)
library(data.table)
library(tidyverse) #提供管道操作
library(sf)
library(raster)
library(terra)
library(spData)
library(ggspatial)

load("~/ClimBehav_bj/results/data/raw/01_AOI/11-北京市-aoi.Rdata")
crs_84 <- st_crs("EPSG:4326")  ## WGS 84

###### Beijing shp ######
beijing <-
  sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/110000_full.json")

#load mp
mp <- fread(file = "~/ClimBehav_bj/results/data/raw/02_MP/MP-beijing-2017.csv.gz", nThread = 14)

timeline <- seq(from = ymd_h(2017010100), to = ymd_h(2017123123), by = "1 hour") %>% 
  as.data.table() %>% 
  transmute(
    time =as.character(year(x) * 1e6 + month(x) * 1e4 + day(x) * 1e2 + hour(x)) 
  )
timeline <- timeline$time

#change data to id - time - cases
checkins.lonlat <- mp %>% 
  pivot_longer(any_of(timeline), names_to = "time", values_to = "visits") %>%
  dplyr::select(id, visits)

beijing_poi <- aoi %>% 
  as_tibble() %>% 
  separate(col = "location", into = c("lon", "lat"), sep = ",", convert = T) %>% 
  dplyr::select(id, lon, lat, type, name) %>%                # type and name for classification
  distinct(id, .keep_all = T)                          #剔除重复的id,提出后为23674条观测

checkins.lonlat <- left_join(checkins.lonlat, beijing_poi, by = "id") %>% as.data.table()

checkins_grid <- 
  # 按照经纬度进行切片
  checkins.lonlat[, ":="(lon_cut = cut(lon,
                                       breaks = seq(115, 118, 0.01),
                                       include.lowest = T),
                         lat_cut = cut(lat,
                                       breaks = seq(39, 41, 0.01),
                                       include.lowest = T))] %>% 
  #统计每个切片的签到数量
  .[, .(checkins_num = sum(visits, na.rm = T)), by = .(lon_cut, lat_cut)] %>% 
  .[, checkins_num_log := log(checkins_num + 1, base = 10)] %>% 
  #计算每个切片的中心点经纬度
  .[, ":="(lon_grid = as.double(lon_cut)*0.01 + 115.005,
           lat_grid = as.double(lat_cut)*0.01 + 39.005)]

checkins_raster <- 
  #生成raster, 使用经纬度和属性值
  rasterFromXYZ(xyz = checkins_grid[, . (lon_grid, lat_grid, checkins_num_log)],
                res = c(0.01, 0.01),
                crs = "EPSG:4326") %>% 
  #mask操作
  raster::mask(beijing) %>% 
  stars::st_as_stars()
#转换为投影坐标系（wgs84）
save(checkins_raster,file = "~/ClimBehav_bj/results/data/processed/checkins_raster.Rdata")
