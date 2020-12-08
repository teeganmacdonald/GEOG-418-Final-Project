install.packages("spgwr")

library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(grid)
library(gridExtra)
library(gtable)
library(sp)


dir <- "C:\\Users\\User\\Documents\\GEOG_418\\FINAL PROJECT"
setwd(dir)
getwd()

elev <- readOGR(dsn = ".","ElevSample") #Read in data
elev <- spTransform(elev, CRS("+init=epsg:26910"))

#Reading in VRI data
VRI <- readOGR(dsn = ".","WatershedVRI") #Read in shapefile
VRI <- spTransform(VRI, CRS("+init=epsg:26910"))
head(VRI@data)

vriCleanCols <- c("FID_VEG_CO", "POLYGON_ID", "PROJ_AGE_1",
                  "SITE_INDEX", "SPECIES__4", "SPECIES__5",
                  "PROJ_HEI_1", "SPECIES_PC", "SPECIES__6",
                  "VRI_LIVE_S", "BASAL_AREA", "WHOLE_STEM",
                  "CROWN_CL_1")

vriClean <- VRI[,vriCleanCols] 
newNames <- c("FID", "PolyID", "Stand_Age", "Site_Index",
              "CoDom_Sp", "Dom_Sp", "Stand_HT", "DomSP_Perc", 
              "CDomSP_Perc", "Stand_Dens", "Stand_BA", "Stand_StemBio", "Stand_CrownCl")

colnames(vriClean@data) <- newNames
vriClean <- vriClean[!is.na(vriClean@data$Stand_StemBio), ]



map_Bio <- tm_shape(vriClean) +
  tm_polygons(col = "Stand_StemBio",
              title = "Stemwood Biomass (tonnes/ha)",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

map_Bio



#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(elev, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(elev)

########### Descriptive Stats ###############
meanElev <- mean(elev$grid_code, na.rm = TRUE)
sdElev <- sd(elev$grid_code, na.rm = TRUE)
medianElev <- median(elev$grid_code, na.rm = TRUE)
modeElev <- as.numeric(names(sort(table(elev$grid_code), decreasing = TRUE))[1])
skewElev <- skewness(elev$grid_code, na.rm = TRUE)
CoVElev <- (sdElev/meanElev)*100
KurtElev <- kurtosis(elev$grid_code, na.rm = TRUE)

meanStem <- mean(vriClean$Stand_StemBio, na.rm = TRUE)
sdStem <- sd(vriClean$Stand_StemBio, na.rm = TRUE)
medianStem <- median(vriClean$Stand_StemBio, na.rm = TRUE)
modeStem <- as.numeric(names(sort(table(vriClean$Stand_StemBio), decreasing = TRUE))[12])
skewStem <- skewness(vriClean$Stand_StemBio, na.rm = TRUE)
CoVStem <- (sdStem/meanStem)*100
KurtStem <- kurtosis(vriClean$Stand_StemBio, na.rm = TRUE)

Samples = c("Elevation", "Stem Biomass") 
Means = c(meanElev, meanStem) 
SD = c(sdElev, sdStem) 
Median = c(medianElev, medianStem) 
Mode <- c(modeElev, modeStem) 
Skewness <- c(skewElev, skewStem) 
CoV <- c(CoVElev, CoVStem)
Kurtosis <- c(KurtElev, KurtStem)

Means <- round(Means, 2)
SD <- round(SD, 2)
Skewness <- round(Skewness, 2)
CoV <- round(CoV, 2)
Kurtosis <- round(Kurtosis, 2)
Median <- round(Median, 2)
Mode <- round(Mode, 2)

data.for.table.1 = data.frame(Samples, Means, SD, Median, Mode)
data.for.table.2 = data.frame(Samples, CoV, Skewness, Kurtosis)

table1 <- tableGrob(data.for.table.1, rows = c("",""))  
t1Caption <- textGrob("Table 1: Central Tendencies for Elevation and Stem Biomass Data", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table.1) + 1)

grid.arrange(table1, newpage = TRUE)
png("Output_Table1.png") 
grid.arrange(table1, newpage = TRUE)
dev.off()

table2 <- tableGrob(data.for.table.2, rows = c("",""))
t2Caption <- textGrob("Table 2: Descriptive Stats for Elevation and Stem Biomass Data", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table2 <- gtable_add_rows(table2, 
                          heights = grobHeight(t2Caption) + padding, 
                          pos = 0)

table2 <- gtable_add_grob(table2,
                          t2Caption, t = 1, l = 2, r = ncol(data.for.table.2) + 1)

grid.arrange(table2, newpage = TRUE)
png("Output_Table2.png") 
grid.arrange(table2, newpage = TRUE)
dev.off()
########## Spatial Autocorrelation ####################

vri.nb <- poly2nb(vriClean)
vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
mi <- moran.test(vriClean$Stand_StemBio, vri.lw, zero.policy = TRUE)
mi


moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(vri.lw)


mI <- mi$estimate[[1]] 
eI <- mi$estimate[[2]] 
var <- mi$estimate[[3]] 

z <- ((mI - eI)/(sqrt(var)))

lisa.test <- localmoran(vriClean$Stand_StemBio, vri.lw, zero.policy = TRUE)

vriClean$Ii <- lisa.test[,1]
vriClean$E.Ii<- lisa.test[,2]
vriClean$Var.Ii<- lisa.test[,3]
vriClean$Z.Ii<- lisa.test[,4]
vriClean$P<- lisa.test[,5]



map_LISA <- tm_shape(vriClean) + 
  tm_polygons(col = "Z.Ii", 
              title = "Local Moran's I", 
              style = "fixed", 
              palette = "viridis", n=6, midpoint = NA, breaks = c(-11, -1.96, 1.96, 37 )) 


map_LISA 

moran.plot(vriClean$Stand_StemBio, vri.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Biomass", 
           ylab="Neighbouring Biomass", quiet=NULL, main = "Local Moran's I Plot for Stem Biomass")

####################### Interpolation ###################################

grd <- as.data.frame(spsample(elev, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(elev)

f.0 <- as.formula(grid_code ~ 1)
var.smpl <- variogram(f.0, elev, cloud = FALSE) 
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill= 37000, model="Exp", range = 4000, nugget=0))
plot(var.smpl, dat.fit, main = "Semivariogram for Elevation")


dat.krg <- krige(f.0, elev, grd, dat.fit)

r <- raster(dat.krg)
r.m <- mask(r, vriClean)

tm_shape(r.m) + 
  tm_raster(n=10, palette="-RdBu",  
            title="Predicted elevation \n(in m)") +
  tm_shape(elev) + tm_dots(size=0) +
  tm_legend(legend.outside=TRUE)+
  tm_layout(title = "Ordinary Kriging Interpolation")


################Combining Data######################

r <- raster(dat.krg)
sufaceMap <- tm_shape(r) + 
  tm_raster(n=5,palette = "viridis",
            title="Elev (m)") +
  tm_shape(elev) + tm_dots(size=0.2)
sufaceMap

vriClean$Elev <- extract(r, vriClean, fun = mean)[,1]

###################Regression##########################


vriClean.no0 <-  vriClean[which(vriClean$Stand_StemBio > 0), ]
vriClean.no0 <-  vriClean.no0[which(vriClean.no0$Elev > 0), ]

plot(vriClean.no0$Stand_StemBio ~ vriClean.no0$Elev)

lm.model <- lm(vriClean.no0$Stand_StemBio ~ vriClean.no0$Elev)
plot(vriClean.no0$Stand_StemBio ~ vriClean.no0$Elev)
abline(lm.model, col = "red")

summary(lm.model)

vriClean.no0$predictlm <- lm.model$fitted.values
vriClean.no0$residuals <- residuals.lm(lm.model)

head(vriClean.no0@data)

map_resid <- tm_shape(vriClean.no0) +
  tm_polygons(col = "residuals",
              title = "Stand Biomass Residuals",
              style = "jenks",
              palette = "viridis", n = 6)
map_resid

res.nb <- poly2nb(vriClean.no0)
res.lw <- nb2listw(res.nb, zero.policy = TRUE, style = "W")
mi <- moran.test(vriClean.no0$residuals, res.lw, zero.policy = TRUE)
mi

mI <- mi$estimate[[1]] 
eI <- mi$estimate[[2]] 
var <- mi$estimate[[3]] 


z <- ((mI - eI)/(sqrt(var)))



######################GWR###########################

vriClean.no0.coords <- sp::coordinates(vriClean.no0)
head(vriClean.no0.coords)

vriClean.no0$X <- vriClean.no0.coords[,1]
vriClean.no0$Y <- vriClean.no0.coords[,2]

GWRbandwidth <- gwr.sel(vriClean.no0$Stand_StemBio ~ vriClean.no0$Elev, 
                        data = vriClean.no0, coords=cbind(vriClean.no0$X,vriClean.no0$Y),adapt=T) 

gwr.model = gwr(vriClean.no0$Stand_StemBio ~ vriClean.no0$Elev, 
                data=vriClean.no0, coords=cbind(vriClean.no0$X,vriClean.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 


gwr.model

results<-as.data.frame(gwr.model$SDF)
head(results)

vriClean.no0$localr <- results$localR2

map_r2 <- tm_shape(vriClean.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style = "jenks",
              palette = "viridis", n = 6)
map_r2

vriClean.no0$coeff <- results$vriClean.no0.Elev
map_coef <- tm_shape(vriClean.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              palette = "viridis", n = 6)
map_coef

################Elevation Distribution##################

kma <- vriClean.no0$Elev
kma$x <- vriClean.no0.coords[,1]
kma$y <- vriClean.no0.coords[,2]


kma.ext <- as.matrix(extent(kma)) 
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)
plot(kma.ppp)

nearestNeighbour <- (nndist(kma.ppp))/1000
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))    
colnames(nearestNeighbour) = "Distance"  
nnd = ((sum(nearestNeighbour$Distance))/N)
nnd <- round(nnd, 3)
elev$area_sqkm <- 205.5
studyArea <- elev$area_sqkm[1]

pointDensity <- (N)/(studyArea)

sqrt(pointDensity)
r.nnd = (1/(2)*(sqrt(pointDensity)))

d.nnd = 1.07453/(sqrt(pointDensity))

R = nnd/r.nnd

SE.NND <- (0.26136/(sqrt(N*pointDensity)))

z = ((nnd - r.nnd)/(SE.NND))


