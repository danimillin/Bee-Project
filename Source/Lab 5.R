#set wd 
setwd("/Users/israel/Desktop/OneDrive - Lawrence University/Biogeography 420/Labs/Lab 5")
library (dismo); library(raster);library(rgdal);library (sp);
library (maptools); library (mgcv); library (gbm)

xy<-read.csv ('V_lagopes_xy.csv')
xy<-xy[,c(3,2)]
duplicated (xy)

#lets remove any duplicate entries why should we do this? 
xy<-xy[!duplicated(xy), ]
duplicated (xy)
str (xy)

#load a simple world map from maptools 
data ("wrld_simpl")

plot (wrld_simpl, xlim=c(-180, -50), ylim=c(23,90), axes=T, col="lightgrey")
points (data=xy, x=xy$longitude, y=xy$latitude, col="red", pch=16)

#creating a background dataset: This is NOT an absence dataset, it simply characterisez the ambient conditionso of the species extent 

#read in our predictor rasters and stack em' 
temp<-raster ("bio1.bil")
temp_cm<- raster ("bio6.bil")
precip<- raster ("bio12.bil")
precip_wq<- raster ("bio13.bil")

predictors<-stack ( temp, temp_cm, precip, precip_wq)

#crop the stack 
predictors_crop<- crop (predictors,  extent(-180,-50,23,80))
plot (predictors_crop)

# we use the first file to create a RasterLayer
plot(predictors_crop[[1]]) 
mask <- raster(predictors_crop[[1]])
mask

# select 500 random points
# set seed to assure that the examples will always have the same random sample.
set.seed(2000)
bg <- randomPoints(mask, 500 )
head (bg)
str (bg)

# set up the plotting area
plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)
# now we repeat the sampling, but limit
# the area of sampling using a spatial extent
e <- extent(-180, -50, 23, 80)
bg2 <- randomPoints(mask, 50, ext=e)
plot(!is.na(mask), legend=FALSE)
plot(e, add=TRUE, col="red")
points(bg2, cex=0.5)

plot (wrld_simpl, xlim=c(-180, -50), ylim=c(23,90), axes=T, col="lightgrey")
points (data=xy, x=xy$longitude, y=xy$latitude, col="red", pch=16, cex=.5)
points (bg, pch=16, col="blue", cex=.25)

#We can now use the extract function as we did before to mine data from each of these geographic locations: 
#presences
presvals <- extract(predictors_crop, xy)
presvals <- na.omit (presvals) #gets rid of any observations with no data 
#psuedo-absences
absvals <- extract(predictors_crop, bg)
absvals <- na.omit (absvals) #gets rid of any observations with no data 

#combine presences and absences into a single dataframe 
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
head(sdmdata)

#investigate co-linearity: Are these predictors telling us the same thing? 
pairs (sdmdata [,2:5], cex=.25, fig=TRUE)

#When two variables are strongly correlated it is best to run the model with only one of them, unless you have a strong biological reason to not do this 

#now our dataset is readly lets fit our first model (Chapter 5 in SDM)
#lets start with a simple Generalized Linear Model (glm)
m1 <- glm(pb ~ bio1 + bio12 + bio13, data=sdmdata)
summary (m1)
response(m1)

#with our model built we can now predict to the entire study extent: 
glm_pred<- predict (predictors_crop, m1)
plot (glm_pred, col=rev(rainbow(50)), main="GLM")# This now plots the habitat suitability for this species 

#evaluate the model 
mE1 <- evaluate(xy, bg, x= predictors_crop , m1)
mE1

#Apply a threshold to convert to a binary map:
?threshold
tr <- threshold(mE1,'spec_sens')
tr
current_pa<- (glm_pred > tr)
plot(glm_pred > tr, main='presence/absence')
#how does this map change if we change the threshold type from spec_sens to the other thresholds? 

#lets now construct a generalized additive model (gam)
library (mgcv)
m2<-gam(pb ~ bio1 * bio12 * bio13, data=sdmdata)
summary (m2)
response(m2)
gam_pred<- predict (predictors_crop, m2)
plot (gam_pred, col=rev(rainbow(50)), main="GAM") 


#Begin Lab 5 here: 
#project our models to a future climate 
#load in rasters and stack them: 
#Important to make sure the files in your future stack are the named the same as your predictor variables 
bio1f<-raster ("bio1.tif")
bio12f<- raster("bio12.tif")
bio13f<- raster("bio13.tif")

future_stack<-stack (bio1f, bio12f, bio13f)
plot (future_stack)

#lets clip it to our study extent: 
future_crop<- crop (future_stack,  extent(-180,-50,23,80))
plot (future_crop)

#now all we do is project our models (m1 and m2) to the future stack! 
summary (m1)
summary (m2)

glm_future<- predict (future_crop, m1)

#you can also apply teh threshold as you before to create a presence absence map 
glm_future_bi<-glm_future > tr
plot (glm_future_bi)

#now plot them

par(mfrow=c(2,2)) 
plot (glm_pred, main="2017 Distribution")
plot (glm_future, main="2070 RCP 8.5")

#using raster math we can calculate the differences between the two maps 
diff_map<- glm_future-glm_pred
plot (diff_map, col=heat.colors(100), main="Difference Map")

diff_pa_map<- glm_future_bi - current_pa
plot (diff_pa_map, col=heat.colors(100), main="Difference Map P/A")

#How would these patterns change if you changed the input scenarios to less dramatic options? 
#BOOM SHAKA-LAKA 

#Lets try a more advanced modeling algorithm: Boosted Regression Trees
head (sdmdata)
summary (sdmdata)
par(mfrow=c(1,1)) 
vulpes_brt<- gbm.step (data=sdmdata, gbm.x= 2:5, gbm.y = 1,
                       family = "bernoulli", tree.complexity = 5,                        
                        learning.rate = 0.01, bag.fraction = 0.5)
#show us the importance of each variable in the model 
summary (vulpes_brt)  

#Plot out the response curves 
gbm.plot (vulpes_brt, n.plots= 4)
  
#compare these fitted curves with those of the glm: Whats the differences? 
response(m1)

#evaluate the model 
pres <- sdmdata[sdmdata[,1]==1, 1]
abs <- sdmdata[sdmdata[,1]==0, 1]
e <- evaluate(p=pres, a=abs)
e

#now use the same predict functions to get the current ranges and future ranges 
brt_pred<- predict(predictors_crop,vulpes_brt,  
                    n.trees=vulpes_brt$gbm.call$best.trees, type="response")
plot (brt_pred)

#evaluate a threshold 
brt_tr = 0.5
plot(brt_pred >brt_tr, main='BRT presence/absence')

bio1f<-raster ("Future rasters/bio1.tif")
bio6f<-raster ("Future rasters/bio6.tif")
bio12f<- raster("Future rasters/bio12.tif")
bio13f<- raster("Future rasters/bio13.tif")

future_stack<-stack (bio1f, bio6f, bio12f, bio13f)

#lets clip it to our study extent: 
future_crop<- crop (future_stack,  extent(-180,-50,23,80))

#now project to a stack of future conditions 
brt_future<- predict(future_crop,vulpes_brt,  
                   n.trees=vulpes_brt$gbm.call$best.trees, type="response")

plot (brt_future, main="2070 BRT RCP 8.5")

par(mfrow=c(2,2)) 
plot (brt_pred, main="2017 Distribution")
plot (brt_future, main="2070 BRT RCP 8.5")

#using raster math we can calculate the differences between the two maps 
diff_map_brt<- brt_future-brt_pred
plot (diff_map_brt, col=heat.colors(100), main="Difference Map")

#Fill in the blank how do we project the binary presence absence maps! 
diff_map_brt_pa<- (brt_future > brt_tr) - (brt_pred > brt_tr) 
plot (diff_map_brt_pa)
