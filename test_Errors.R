library(yamdar)
library(superheat)


load("errordat.Rdata")

modMat = lapply(modH, toHypotMatrix)
x = Yamda(data = dat, modMat, factors = T, nneg = T)
x$stats

plot_expected(original_matrix = cor(dat), yamda_object = x, n=1)

(x = Yamda(data = dat, modMat, factors = T, nneg = F))$stats
(x = Yamda(data = dat, modMat, factors = F, nneg = T))$stats
(x = Yamda(data = dat, modMat, factors = T, nneg = F))$stats
(x = Yamda(data = dat, modMat, factors = F, nneg = F))$stats

tryCatch(expr = {Yamda(data = dat, modMat, factors = F, nneg = F)}, 
             error = function(e) return(NA))

calcZTransCoef(modMat[[1]], cor(dat), F)


data = dat
hypot = modMat[[4]]
nneg = F
factors = F

fitModuleCoef(data = dat, hypot =  modMat[[1]], nneg = T, factors = F)
fitModuleCoef(data = dat, hypot =  modMat[[1]], nneg = T, factors = T)
fitModuleCoef(data = dat, hypot =  modMat[[1]], nneg = F, factors = F)
fitModuleCoef(data = dat, hypot =  modMat[[1]], nneg = F, factors = T)

load("./exemplo.Rdata")

h_list
out = Yamda(simdata[[9]], h_list)
out$stat

Yamda(simdata[[9]],h_list)$stat
Yamda(simdata[[9]],h_list[1])$stat

plot_expected(cor(simdata[[9]]), out, n=2, ncol = 3)


#First Matrix
b<-rnorm(16,0.6,0.01) #between modules covariances
# m<-rnorm(45,0.65,0.01) #within modules covariances
m<-rnorm(32,0.60,0.01) #within modules covariances/ No modularity
M1<-diag(16)
M1[9:16,1:8]<-b
M1[lower.tri(M1) & !M1>0]<-m
M1[upper.tri(M1)]<-t(M1)[upper.tri(M1)]
M1<-nearPD(M1)$mat %>% as.matrix()
rmvnorm(45,sigma = M1)