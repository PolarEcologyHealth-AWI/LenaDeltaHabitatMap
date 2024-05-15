### Hyperspectral
library(tidyverse)
library(vegan)

### Reflectance / S2 Bands
trans <- tibble(S2 = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12"),
                lower = c(433, 458, 543, 650, 698, 733, 773, 785, 855, 935, 1360, 1565, 2100),
                upper = c(453, 523, 578, 680, 713, 748, 793, 900, 875, 955, 1390, 1655, 2280))


fls <- list.files("Data/Reflectance-Spectra/", pattern = "Final", full.names = T)

tab <- do.call("rbind", parallel::mclapply(fls, function(f) {
  out <- suppressMessages(read_csv(f, progress = F)) %>% 
    rename(id = sub_plot) %>% dplyr::select(wavelength, id, meas, percent_reflectance) %>%
    mutate(S2 = unlist(sapply(wavelength, function(x) {
      out <- trans$S2[x>=trans$lower & x<=trans$upper]
      if(length(out)>0) out[1] else NA
    }))) %>% filter(!is.na(S2))
  
  out %>% group_split(meas) %>% lapply(., function(s) {
    s %>% dplyr::select(meas, S2, percent_reflectance) %>% 
      pivot_wider(names_from = S2, values_from = percent_reflectance, values_fn = median) %>%
      mutate(ndvi = (B8-B4)/(B8+B4),
             id   = s$id[1]) %>% relocate(id, .before = meas) %>% dplyr::select(-B10)
  }) %>% do.call("rbind", .)
}, mc.cores = 9))

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
                          '#d2042d', '#f19414', '#202cff'),
                  sim = c(10, 2, 1, 11, 7, 5, 9, 6, 4, 3, 8, 12, 13)-1)

classesPlots <- tibble(id    = 1:27,
                       class = c(5, 5, 8, 1, 6, 6, 5, 9, 5, 5, 6, NA, 8, 5, 4, 1, 0, 9, NA, 6, 4, 5, 4, NA, NA, 5, 9)) %>%
  left_join(cls, by = "id")

tabClasses <- tab %>% left_join(classesPlots %>% dplyr::select(id, class), by = "id") %>% 
  filter(!is.na(class)) %>% relocate(class, .before = meas)

veg.mds     <- metaMDS(tabClasses[,-c(1:3)], distance = "bray", autotransform = FALSE)
veg.spp.fit <- envfit(veg.mds, tabClasses[,-c(1:3)], permutations = 999)

plot(NA, xlim = c(-0.45, 0.3), ylim = c(-0.3, 0.2), type = "n", las = 1, xlab = "NMDS 1 []", ylab = "NMDS 2 []")
orditorp(veg.mds, display = "sites", labels = F, 
         pch = 16, col = adjustcolor("grey50", alpha.f = 0.4), cex = 1)
ordiellipse(veg.mds, groups = tabClasses %>% pull(class), draw = "polygon", lty = 1, 
            col = cls$cls[match(as.numeric(names(table(tabClasses[,2]))), cls$id)], alpha = 0.7)
# plot(veg.mds)
plot(veg.spp.fit, p.max = 0.01, col = "black", cex = 0.7)
legend("bottomleft", pch = 22, pt.cex = 3, pt.bg = cls$cls[1:11], cls$names[1:11], cex = 1.2, bg = adjustcolor("white", alpha.f = 0.5), box.col = NA)
# dev.off()

dendroTab <- tabClasses[,-c(1,3)] %>% group_by(class) %>% summarise_all(median) %>% as.data.frame()
rownames(dendroTab) <- cls$names[match(dendroTab$class, cls$id)]
hc <- hclust(dist(dendroTab))
plot(hc, hang = -1)




### Overlap

cls <- cls %>% arrange(sim)

# scores
spec   <- c(0:10)
scNMDS <- tibble(class = tabClasses$class, scores(veg.mds)$sites %>% as.data.frame())
D_mat  <- matrix(nrow = length(0:10), ncol = length(0:10))

ind    <- data.frame(expand.grid(0:10, 0:10), true = c(lower.tri(D_mat)))
D_out  <- abind::abind(D_mat, D_mat, D_mat, along = 3)

for(s in which(ind$true)) {
  
  class1 <- cls$id[cls$sim==ind[s,1]]
  class2 <- cls$id[cls$sim==ind[s,2]]
  
  extrT <- as.data.frame(scNMDS) %>% as_tibble() %>%
    filter(class%in%c(class1, class2))
  
  if(class1%in%extrT$class & class2%in%extrT$class) {
  ## pca
  pca.cal <- ade4::dudi.pca(extrT[,-1], center = T, scale = T, scannf = F, nf = 2)
  
  ### Grid
  z1<- grid.clim(pca.cal$li, pca.cal$li[extrT$class==class1,], pca.cal$li[extrT$class==class1,], R = 100)
  z2<- grid.clim(pca.cal$li, pca.cal$li[extrT$class==class2,], pca.cal$li[extrT$class==class2,], R = 100)
  
  ### Niche overlap/equivalency test
  # eT <- equivalency.test(z1, z2, 100, 10)
  eT <- c(niche.overlap(z1, z2, cor=T)$D, NA)
  # sT1 <- similarity.test(z1, z2,  100, 4)
  # sT2 <- similarity.test(z2, z1,  100, 4)
  
  D_out[which(spec==ind[s,1]), which(spec==ind[s,2]),] <- c(eT[1], eT[2], NA)
  }
}
# save(D_out, file = "D_out.rda")
# load("D_out.rda")


pols1 <- rasterToPolygons(raster(t(D_out[,,1])))
plot(pols1, col = rev(viridis::mako(11))[cut(pols1$layer*100, seq(0, 100, 10), labels = FALSE)], border = NA)
tmp <- tibble(coordinates(pols1) %>% as.data.frame(), dat = pols1$layer) %>% filter(!is.na(dat))
text(tmp$V1, tmp$V2, round(tmp$dat*100,0), cex = 2)
axis(1)
