#Code for surface mass balance prediction as in Section 3 of "A review of Bayesian modelling in glaciology"
#Uses R-INLA to make spatial predictions and infer uncertainties of surface mass balance

#If you don't have the following R packages, you can install with install.packages
#However, INLA is not on CRAN but can be downloaded with 
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#Check https://www.r-inla.org/ for most recent download and installtion instructions.

library(INLA)
library(splancs)
library(plyr)

#NOTE ON ASSUMED VARIABLES for R script
# - PRborder: a 2 column matrix of x and y coordinates of border surrounding region to be modelled. (For instance, rectangular.)
# - x_vals: vector of x coordinates in the grid being modeled
# - y_vals: vector of y coordinates in the grid being modeled
# - coords: 2 column matrix of x and y coordinates of the locations of the observed surface mass balance data
# - Y: vector that contains observed mass balance measurement values for each of the coordinates in coords
# - elevation: vector that contains surface elevation values for each of the coordinates in coords
# - elevation.prd: vector that contains the surface elevation values for the to-be-predicted x, y coordinates
#Suggestion: scale x and y coordinates to help R-INLA calculations go faster

#create a mesh for SPDE method -- coords contains an n by 2 matrix for x (column) and y coordinates 
prmesh <- inla.mesh.2d(coords,loc.domain=PRborder,max.edge=.1)
nxy <- c(length(x_vals),length(y_vals))
projgrid <- inla.mesh.projector(prmesh,xlim=range(PRborder[,1]),ylim=range(PRborder[,2]),dims=nxy)
xy.in <- inout(projgrid$lattice$loc,cbind(PRborder[,1],PRborder[,2]))

  
#observation matrix A
A <- inla.spde.make.A(prmesh,loc=coords)
#Note: prior.range and prior.sigma specify PC prior hyperparameters
spde <- inla.spde2.pcmatern(prmesh,alpha=2,prior.range = c(.002,.01),prior.sigma = c(.3,.05))


#inla.stack function
mesh.index <- inla.spde.make.index(name='field',n.spde=spde$n.spde)
stk.dat <- inla.stack(data=list(y=Y),A=list(A,1),tag="est",
                        effects=list(c(mesh.index,list(Intercept=1)),
                                     list(long=inla.group(coords[,1]),
                                          lat=inla.group(coords[,2]),
                                          elevation=inla.group(elevation))))

#model formula
f.s <- y ~ -1 + Intercept + long + lat + elevation + f(field, model = spde)
coord.prd <- projgrid$lattice$loc[xy.in,]
A.prd <- projgrid$proj$A[xy.in,]
ef.prd = list(c(mesh.index,list(Intercept=1)),list(long=inla.group(coord.prd[,1]),lat=inla.group(coord.prd[,2]),elevation=inla.group(elevation.prd)))
stk.prd <- inla.stack(data=list(y=NA),A=list(A.prd,1),tag="prd",effects=ef.prd)
stk.all <- inla.stack(stk.dat,stk.prd)

r2.s <- inla(f.s,family='Normal',
               data=inla.stack.data(stk.all),
               control.fixed=list(prec=list(elevation=.1,long=1,lat=1)),
               control.predictor=list(A=inla.stack.A(stk.all),
                                      compute=TRUE,link=1),
               quantiles=NULL,
               verbose=TRUE,
               control.inla = list(int.strategy = "eb"))

#extract indices for prediction nodes and then the mean and standard deviation of the response
id.prd <- inla.stack.index(stk.all,"prd")$data
sd.prd <- m.prd <- matrix(NA,nxy[1],nxy[2])
#mean prediction of mass balance over grid points
m.prd[xy.in] <- r2.s$summary.fitted.values$mean[id.prd]
#standard deviation of mass balance prediction over grid points
sd.prd[xy.in] <- r2.s$summary.fitted.values$sd[id.prd]
