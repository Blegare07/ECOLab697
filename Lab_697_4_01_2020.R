#Mini_Lab
if(!require("rspatial")) devtools::install_github('rspatial/rspatial')

library(rspatial)
d <- sp_data('precipitation')
head(d)

class(d)


d$prec <- rowSums(d[, c(6:17)])

plot(
  sort(d$prec),
  main = "CA annual precipitation",
  ylab = 'Annual precipitation (mm)', xlab = 'Stations',
  las = 1 )


require(sp)
dsp <- SpatialPoints(d[,4:3], proj4string = CRS("+proj=longlat +datum=NAD83"))
dsp <- SpatialPointsDataFrame(dsp, d)
CA <- sp_data("counties")

cuts<-c(0,200,300,500,1000,3000)

blues <- colorRampPalette(c('yellow', 'orange', 'blue', 'dark blue'))
pols <- list("sp.polygons", CA, fill = "lightgray")
spplot(dsp, 'prec', cuts = cuts, col.regions = blues(5), sp.layout = pols, pch = 20, cex = 2)

TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0
+lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83
+units=m +ellps=GRS80 +towgs84=0,0,0")
require(rgdal)
dta <- spTransform(dsp, TA)
cata <- spTransform(CA, TA)


RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm = TRUE))
}


null <- RMSE(mean(dsp$prec), dsp$prec)
null



#Proximity Polygon Interpolation

require(dismo)


v <- voronoi(dta)
plot(v)

require(rgeos)
ca <- aggregate(cata)
vca <- intersect(v, ca)
spplot(vca, 'prec', col.regions = rev(get_col_regions()))

ca
cata


par(mfrow = c(1, 2)); plot(cata, main = "cata"); plot(ca, main = "ca")


r <- raster(cata, res = 10000)
vr <- rasterize(vca, r, 'prec')
plot(vr)


set.seed(5132015)
kf <- kfold(nrow(dta))
rmse <- rep(NA, 5)
for (k in 1:5) {
  test <- dta[kf == k, ]
  train <- dta[kf != k, ]
  v <- voronoi(train)
  p <- extract(v, test)
  rmse[k] <- RMSE(test$prec, p$prec)
}
rmse


mean(rmse)
1-(mean(rmse)/null)

#Nearest neighbour interpolation


require(gstat)



gs <- gstat(
  formula = prec~1,
  locations = dta,
  nmax = 5,
  set = list(idp = 0))
nn <- interpolate(r, gs)


nnmsk <- mask(nn, vr)
plot(nnmsk)



rmsenn <- rep(NA, 5)
for (k in 1:5) {
  test <- dta[kf == k, ]
  train <- dta[kf != k, ]
  gscv <- gstat(
    formula = prec~1,
    locations = train,
    nmax = 5,
    set = list(idp = 0))
  p <- predict(gscv, test)$var1.pred
  rmsenn[k] <- RMSE(test$prec, p)
}

rmsenn


mean(rmsenn)


1 - (mean(rmsenn) / null)


require(gstat)
gs <- gstat(formula = prec~1, locations = dta)
idw <- interpolate(r, gs)


idwr <- mask(idw, vr)
plot(idwr)

rmse <- rep(NA, 5)
for (k in 1:5) {
  test <- dta[kf == k, ]
  train <- dta[kf != k, ]
  gs <- gstat(formula = prec~1, locations = train)
  p <- predict(gs, test)
  rmse[k] <- RMSE(test$prec, p$var1.pred)
}
rmse
mean(rmse)
1 - (mean(rmse) / null)


gs2 <- gstat(formula = prec~1, locations = dta, nmax = 1, set = list(idp = 1))


require(rspatial)
x <- sp_data("airqual")
x$OZDLYAV <- x$OZDLYAV * 1000



require(sp)
require(rgdal)
coordinates(x) <- ~LONGITUDE + LATITUDE
proj4string(x) <- CRS('+proj=longlat +datum=NAD83')
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +ellps=GRS80")
aq <- spTransform(x, TA)


cageo <- sp_data('counties.rds')
ca <- spTransform(cageo, TA)
r <- raster(ca)
res(r) <- 10 # 10 km if your CRS's units are in km
g <- as(r, 'SpatialGrid')


#fit a variogram

require(gstat)
gs <- gstat(formula = OZDLYAV~1, locations = aq)
v <- variogram(gs, width = 20)
head(v)


plot(v)


fve <- fit.variogram(v, vgm(85, "Exp", 75, 20))
fve


plot(variogramLine(fve, 400), type = 'l', ylim = c(0,120))
points(v[,2:3], pch = 20, col = 'red')


fvs <- fit.variogram(v, vgm(85, "Sph", 75, 20))
fvs


plot(variogramLine(fvs, 400), type = 'l', ylim = c(0,120) ,col = 'blue', lwd = 2)
points(v[,2:3], pch = 20, col = 'red')

plot(v, fve)

#ordinary kriging

k <- gstat(formula = OZDLYAV~1, locations = aq, model = fve)
# predicted values
kp <- predict(k, g)

spplot(kp)


ok <- brick(kp)
ok <- mask(ok, ca)
names(ok) <- c('prediction', 'variance')
plot(ok)


#compare with other methods

require(gstat)
idm <- gstat(formula = OZDLYAV~1, locations = aq)
idp <- interpolate(r, idm)

idp<-mask(idp,ca)
plot(idp)


RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm = TRUE))
}
f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula = OZDLYAV~1, locations = train, nmax = nmx, set = list(idp = idp))
  p <- predict(m, newdata = test, debug.level = 0)$var1.pred
  RMSE(test$OZDLYAV, p)
}
set.seed(20150518)
i <- sample(nrow(aq), 0.2 * nrow(aq))
tst <- aq[i,]
trn <- aq[-i,]
opt <- optim(c(8, .5), f1, test = tst, train = trn)
opt


m <- gstat(formula = OZDLYAV~1, locations = aq, nmax = opt$par[1], set = list(idp = opt$par[2]))
idw <- interpolate(r, m)


idw <- mask(idw, ca)
plot(idw)

#a thin plate spline model

require(fields)

m <- Tps(coordinates(aq), aq$OZDLYAV)
tps <- interpolate(r, m)
tps <- mask(tps, idw)
plot(tps)

#Cross-validate

require(dismo)
nfolds <- 5
k <- kfold(aq, nfolds)
ensrmse <- tpsrmse <- krigrmse <- idwrmse <- rep(NA, 5)
for (i in 1:nfolds)
{
  test <- aq[k != i,]
  train <- aq[k == i,]
  m <- gstat(formula = OZDLYAV~1, locations = train, nmax = opt$par[1], set = list(idp = opt$par[2]))
  p1 <- predict(m, newdata = test, debug.level = 0)$var1.pred
  idwrmse[i] <- RMSE(test$OZDLYAV, p1)
  m <- gstat(formula = OZDLYAV~1, locations = train, model = fve)
  p2 <- predict(m, newdata = test, debug.level = 0)$var1.pred
  krigrmse[i] <- RMSE(test$OZDLYAV, p2)
  m <- Tps(coordinates(train), train$OZDLYAV)
  p3 <- predict(m, coordinates(test))
  tpsrmse[i] <- RMSE(test$OZDLYAV, p3)
  w <- c(idwrmse[i], krigrmse[i], tpsrmse[i])
  weights <- w / sum(w)
  ensemble <- p1 * weights[1] + p2 * weights[2] + p3 * weights[3]
  ensrmse[i] <- RMSE(test$OZDLYAV, ensemble)
}
rmi <- mean(idwrmse)
rmk <- mean(krigrmse)
rmt <- mean(tpsrmse)
rms <- c(rmi, rmt, rmk)
rms

rme <- mean(ensrmse)
rme

weights <- ( rms / sum(rms) )
s <- stack(idw, ok[[1]], tps)
ensemble <- sum(s * weights)


s <- stack(idw, ok[[1]], tps, ensemble)
names(s) <- c('IDW', 'OK', 'TPS', 'Ensemble')
plot(s)
