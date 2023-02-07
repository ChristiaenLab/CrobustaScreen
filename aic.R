source('modelfns.R')
library(keras)

dat <- read.csv('out/z_dat.csv',row.names=1)
dat <- dat[,apply(abs(dat),2,max)!=0]
x_train <- t(t(dat)/apply(abs(dat),2,max))

models <- c("2022-09-28/encode2.model",
	    "2022-09-28/encode3.model",
	    "2022-10-02/encode7.model",
	    "2022-10-03/encode14.model")

models <- lapply(models,
		 load_model_hdf5)
err <- sapply(models,evaluate,x_train,x_train)

layers <- sapply(models,function(x) length(x$layers))

encoded <- mapply(function(x,y) model.out(x,x_train,y/2),models,layers)

aic <- mapply(function(k,L) 2*k-2*log(L),sapply(encoded,ncol),err)

tab <- data.frame(bottleneck=sapply(encoded,ncol),layers=layers,MSE=err,AIC=aic)

write.csv(tab,'aic.csv')
