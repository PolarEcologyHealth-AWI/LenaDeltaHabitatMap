library(tidyverse)
library(stars)
library(sf)
sf_use_s2(FALSE)

# tif <- read_stars("Data/Landgraf_Pangea/Sentinel2_2018_CentralLenaDelta_landcover_Landgraf.tif") %>% st_as_stars() %>%
#   setNames("Classes")
# tif$Classes <- as.factor(tif$Classes+1)
# plot(tif)
# 
# 
# cls <- data.frame(id = 0:12,
#                   names = c('Moist Equisetum and Shrubs on Floodplain',
#                             'Dry Low Shrub Community',
#                             'Moist to Wet Sedge Complex',
#                             'Wet Sedge Complex',
#                             'PC_50%: Wet Polygon Complex',
#                             'PC_20%: Moist Polygon Complex',
#                             'Dry Grass to Wet Sedge Complex',
#                             'Sparsely Vegetated Areas',
#                             'Dry Tussock Tundra',
#                             'PC_10%: Dry Polygon Complex',
#                             'Dwarsh Shrub-herb Communities',
#                             'Sand',
#                             'Water'),
#                   cls = c('white',
#                           '#68ab5f', '#1c5f2c', '#a3cc51', '#43DF4F', '#ccb879',
#                           '#dcd939', '#b5c58f', '#b8d9eb', '#af963c', '#6c9fb8',
#                           'white', 'white'))
# 
# 
# 
# ggplot() +
#   geom_stars(data = tif, downsample = 2) +
#   scale_fill_manual(values = cls$cls) +
#   coord_equal() + 
#   theme_void()


## colors
cls <- data.frame(id = 0:12,
                  names = c('Moist Equisetum and Shrubs on Floodplain',
                            'Dry Low Shrub Community',
                            'Moist to Wet Sedge Complex',
                            'Wet Sedge Complex',
                            'PC_50%: Wet Polygon Complex',
                            'PC_20%: Moist Polygon Complex',
                            'Dry Grass to Wet Sedge Complex',
                            'Sparsely Vegetated Areas',
                            'Dry Tussock Tundra',
                            'PC_10%: Dry Polygon Complex',
                            'Dwarsh Shrub-herb Communities',
                            'Sand',
                            'Water'),
                  cls = c('#68ab5f', '#1c5f2c', '#a3cc51', '#43DF4F', '#ccb879',
                          '#dcd939', '#b5c58f', '#b8d9eb', '#af963c', '#6c9fb8',
                          '#d2042d', 'transparent', 'transparent'))



### Water mask
fls <- list.files("~/Google Drive/My Drive/Science/ProjectsData/LenaDeltaClassification/", pattern = "Sand", full.names = T)
sand <- read_stars(fls[1])
# plot(sand)

bbox1 <- st_bbox(c(xmin = 126.673, xmax = 126.513, ymax = 72.57212, ymin = 72.65112), crs = st_crs(4326)) %>% st_as_sfc()
bbox2 <- st_bbox(c(xmin = 124.7487, xmax = 124.9699, ymax = 73.53919, ymin = 73.60408), crs = st_crs(4326)) %>% st_as_sfc()
bbox3 <- st_bbox(c(xmin = 124.8347, xmax = 124.9883, ymax = 72.78934, ymin = 72.83981), crs = st_crs(4326)) %>% st_as_sfc()

bboxes <- do.call("rbind", list(bbox1, bbox2, bbox3)) %>% st_as_sfc() %>% st_set_crs(4326)

### Overall Map
pls <- read_sf("~/Google Drive/My Drive/Science/ProjectsData/LenaDeltaClassification/LenaDelta_LandCoverHabitatClasses_ROI/LenaDelta_LandCoverHabitatClasses_ROI.shp")
tif <- read_stars("~/Google Drive/My Drive/Science/ProjectsData/LenaDeltaClassification/LenaDelta_LandCoverHabitatClasses.tif") %>%
  setNames("Classes") %>% st_crop(pls) %>% st_as_stars(downsample = 50)
tif$Classes <- as.factor(tif$Classes+1)

ggplot() +
  geom_stars(data = tif, downsample = 0) + 
  scale_fill_manual(values = cls$cls, na.value = "white") +
  # geom_sf(data = bboxes, mapping = aes(geometry = geometry), fill = NA) +
  theme_void()


## inset 1
tif2.1 <- read_stars("/Volumes/GoogleDrive-116402590020799794814/My Drive/Science/ProjectsData/LenaDeltaClassification/LenaDelta_LandCoverHabitatClasses.tif") %>%
  setNames("class") %>% st_crop(bbox1) %>% st_as_stars()
tif2.2 <- sand %>% st_crop(bbox1) %>% st_as_stars() %>% setNames("sand") %>% st_warp(tif2.1)
tif2   <- c(tif2.1, tif2.2) %>% merge() %>% setNames("classes")

tifRast <- as(tif2, "Raster")

raster::plot(tifRast$sand, breaks = seq(0.025, 1, length = 100), col = rev(viridis::inferno(99)), legend = FALSE)
raster::plot(tifRast$class, breaks = -0.5:12.5, col = cls$cls, legend = FALSE, add = T)



## inset 2
tif2.1 <- read_stars("/Volumes/GoogleDrive-116402590020799794814/My Drive/Science/ProjectsData/LenaDeltaClassification/LenaDelta_LandCoverHabitatClasses.tif") %>%
  setNames("class") %>% st_crop(bbox2) %>% st_as_stars()
tif2.2 <- sand %>% st_crop(bbox2) %>% st_as_stars() %>% setNames("sand") %>% st_warp(tif2.1)
tif2   <- c(tif2.1, tif2.2) %>% merge() %>% setNames("classes")

tifRast <- as(tif2, "Raster")


raster::plot(tifRast$sand, breaks = seq(0.025, 1, length = 100), col = rev(viridis::inferno(99)), legend = FALSE)
raster::plot(tifRast$class, breaks = -0.5:12.5, col = cls$cls, legend = FALSE, add = T)

## inset 3
tif2.1 <- read_stars("/Volumes/GoogleDrive-116402590020799794814/My Drive/Science/ProjectsData/LenaDeltaClassification/LenaDelta_LandCoverHabitatClasses.tif") %>%
  setNames("class") %>% st_crop(bbox3) %>% st_as_stars()
tif2.2 <- sand %>% st_crop(bbox3) %>% st_as_stars() %>% setNames("sand") %>% st_warp(tif2.1)
tif2   <- c(tif2.1, tif2.2) %>% merge() %>% setNames("classes")

tifRast <- as(tif2, "Raster")


raster::plot(tifRast$sand, breaks = seq(0.025, 1, length = 100), col = rev(viridis::inferno(99)), legend = FALSE)
raster::plot(tifRast$class, breaks = -0.5:12.5, col = cls$cls, legend = FALSE, add = T)
