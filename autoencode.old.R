source("imageFns.R")
source('modelfns.R')
library(dirfns)
library(moreComplexHeatmap)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, '--threshold', 
		     action = 'store_true')
parser <- add_option(parser, '--dat', action = 'store',
		     default = 'out/z_dat.csv')
parser <- add_option(parser, '--groups', action = 'store',
		     default = 'out/groups.csv')
opts <- parse_args(parser)


groups <- read.csv(opts$groups)
dat <- read.csv(opts$dat,row.names=1)

if(opts$threshold){
	dat <- cbind(dat,zdf(groups[,c("Threshold.nuclei","Threshold.membrane")]))
}

dat <- dat[,apply(abs(dat),2,max)!=0]
x_train <- t(t(dat)/apply(abs(dat),2,max))
cols <- groups[,c("Condition","Phenotype")]

corHeatmap(dat,append.date=T)

if(opts$threshold){
	nuc.sel <- which(x_train[,115]==0)
	mem.sel <- which(x_train[,116]==0)
	sel <- union(nuc.sel,mem.sel)

	denoise.nuc <- denoise(x_train[-sel,-115],
			       x_train[-sel,115],
			       c(115,58,29,14,8,4,2,1),
			       'denoise.nuc')
	denoise.mem <- denoise(x_train[-sel,-116],
			       x_train[-sel,116],
			       c(115,58,29,14,8,4,2,1),
			       'denoise.mem')
	denoise.thr <- denoise(x_train[-sel,-115:-116], 
			       x_train[-sel,115:116], 
			       c(114,58,29,14,8,4,2), 
			       'denoise.thr')

	x_train[mem.sel[1],
		116] <- predict(denoise.mem,
				x_train[mem.sel[1],
					-116,drop=F])
	x_train[nuc.sel,
		115:116] <- predict(denoise.thr,
				    x_train[nuc.sel,
					    -115:-116,
					    drop=F])
}

dir.csv(x_train,'dat')

hm <- Heatmap(t(x_train),name='standardized z-score',
	     #              cell.h=0.024, cell.w=0.005, 
	     column_split=cols$Condition,
	     show_column_names=F,
	     height = unit(.8,"npc"),
	     width=unit(.7,"npc"),
	     column_title_rot = 90)

dir.pdf('params',height=20,width=20)
draw(hm)
dev.off()

encode2 <- autoencode(x_train,c(58,29,14,8,4,2),
		      cols,'encode2',epochs=10000)
encode3 <- autoencode(x_train,c(58,29,14,7,3),
		      cols,'encode3',epochs=10000)
encode7 <- autoencode(x_train,c(58,29,14,7),
		      cols,'encode7',epochs=10000)
encode14 <- autoencode(x_train,c(58,29,14),
		      cols,'encode14',epochs=10000)

models <- list(encode2,encode3,encode7,encode14)

err <- sapply(models,evaluate,x_train,x_train)

layers <- sapply(models,function(x) length(x$layers))

encoded <- mapply(function(x,y) model.out(x,x_train,y/2),models,layers)

aic <- mapply(function(k,L) 2*k-2*log(L),sapply(encoded,ncol),err)

tab <- data.frame(bottleneck=sapply(encoded,ncol),layers=layers,MSE=err,AIC=aic)

dir.csv(tab,'aic')

encoded <- model.out(encode2,x_train,6)
colnames(encoded) <- c('encoding1','encoding2')

clusts2 <- encoded.leiden(encode2,x_train,6,'encode2')
clusts3 <- encoded.leiden(encode3,x_train,5,'encode3')
clusts7 <- encoded.leiden(encode7,x_train,4,'encode7')
clusts14 <- encoded.leiden(encode14,x_train,3,'encode14')


sil2 <- leiden.sil(clusts2$clusts,'2022-10-03/encode2')
sil2.mean <- sapply(sil2,sapply,function(x) mean(x[,3]))
clusts2$log2error*sil2.mean

sil3 <- leiden.sil(clusts3$clusts,'2022-10-03/encode3')
sil3.mean <- sapply(sil3,sapply,function(x) mean(x[,3]))
clusts3$log2error*sil3.mean

sil7 <- leiden.sil(clusts7$clusts,'2022-10-03/encode7')
sil7.mean <- sapply(sil7,sapply,function(x) mean(x[,3]))
clusts7$log2error*sil7.mean

sil14 <- leiden.sil(clusts14$clusts,'2022-10-03/encode14')
sil14.mean <- sapply(sil14,sapply,function(x) mean(x[,3]))
max(clusts14$log2error*sil14.mean)

nclust <- sapply(clusts2[2:11],sapply,max)

dir.pdf('log2error')
Heatmap(clusts2$log2error,cluster_rows = F,cluster_columns = F,name='-log2(error)')
dev.off()

dir.pdf('nclust')
Heatmap(nclust,cluster_rows = F,cluster_columns = F,name='number of clusters')
dev.off()

dir.pdf('err_nclust')
Heatmap(-log2((2^(-clusts2$log2error))/nclust),cluster_rows = F,cluster_columns = F,name='log2(error*numer of clusters)')
dev.off()
clusts <- apply(clusts,2,as.character)
colnames(clusts) <- paste0('k',as.character(seq(3,21,2)))

err <- mapply(test.knn,k=seq(3,21,2),clust=as.data.frame(clusts),MoreArgs = list(dat=encoded))

lab <- paste('k =',as.character(seq(3,21,2)),' -log2(err) =',as.character(round(err,2)))

plot.clust(encoded,as.data.frame(clusts),'leiden0.5',labs=lab,width=20)

mse.plot(encode7,x_train,x_train,clusts$umap,cols,'2022-09-06/encode7',ncols=2)

clusts17 <- sapply(seq(0.1,1,0.1),function(r) get.clusts(encoded,17,r))
clusts17 <- apply(clusts17,2,as.character)
colnames(clusts17) <- paste0('leiden',as.character(seq(0.1,1,0.1)))

err17 <- mapply(test.knn,clust=as.data.frame(clusts17),MoreArgs = list(dat=encoded,k=17))

lab17 <- paste('res =',as.character(seq(0.1,1,0.1)),' -log2(err) =',as.character(round(err17,2)))

plot.clust(encoded,as.data.frame(clusts17),'k17',labs=lab17,width=20)

# set model
model <- keras_model_sequential()
model %>%
  layer_dense(units = 6, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
  layer_dense(units = 6, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

# view model layers
summary(model)

# compile model
model %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam"
)

# fit model
model %>% fit(
  x = x_train, 
  y = x_train, 
  epochs = 2000,
  verbose = 0
)

# evaluate the performance of the model
mse.ae2 <- evaluate(model, x_train, x_train)
mse.ae2

# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, index=6)$output)
intermediate_output <- predict(intermediate_layer_model, x_train)
write.csv(intermediate_output,'out/bottleneck')

# plot the reduced dat set
aedf3 <- data.frame(node1 = intermediate_output[,1], node2 = intermediate_output[,2], node3 = intermediate_output[,3])
ae_plotly <- plot_ly(aedf3, x = ~node1, y = ~node2, z = ~node3, color = ~ais$sex) %>% add_markers()

pdf('bottleneckcond.pdf')
ggplot(data.frame(bottleneck1 = intermediate_output[,1], bottleneck2 = intermediate_output[,2]), aes(x = bottleneck1, y = bottleneck2, col = groups$Condition)) + geom_point()
dev.off()
pdf('bottleneckpheno.pdf')
ggplot(data.frame(bottleneck1 = intermediate_output[,1], bottleneck2 = intermediate_output[,2]), aes(x = bottleneck1, y = bottleneck2, col = groups$Phenotype)) + geom_point()
dev.off()

encoder <- keras_model_sequential()
encoder %>%
  layer_dense(units = 58, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 29, activation = "tanh") %>%
  layer_dense(units = 15, activation = "tanh") %>%
  #   layer_dense(units = 8, activation = "tanh") %>%
  #   layer_dense(units = 4, activation = "tanh") %>%
  #   layer_dense(units = 2, activation = "tanh") %>%
  #   layer_dense(units = 4, activation = "tanh") %>%
  #   layer_dense(units = 8, activation = "tanh") %>%
  #   layer_dense(units = 15, activation = "tanh") %>%
  layer_dense(units = 29, activation = "tanh") %>%
  layer_dense(units = 58, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

# compile model
encoder %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam"
)

# fit model
encoder %>% fit(
  x = x_train, 
  y = x_train, 
  epochs = 10000,
  verbose = 0
)

encoder.mse <- evaluate(encoder, x_train, x_train)

encoder.out <- keras_model(inputs = encoder$input, outputs = get_layer(encoder, index=3)$output)
encoder.bottleneck <- predict(encoder.out, x_train)
encoder.umap <- umap(encoder.bottleneck)

condplot <- ggplot(as.data.frame(encoder.umap$layout), aes(x=V1,y=V2,col = groups$Condition)) + geom_point()
phenoplot <- ggplot(as.data.frame(encoder.umap$layout), aes(x = V1, y = V2, col = groups$Phenotype)) + geom_point()

g <- ggarrange(condplot,phenoplot,labels=c('condition','phenotype'))

ggsave('encoder.pdf',g,width = 20)

encoder2 <- keras_model_sequential()
encoder2 %>%
  layer_dense(units = 58, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 29, activation = "tanh") %>%
  layer_dense(units = 15, activation = "tanh") %>%
  layer_dense(units = 8, activation = "tanh") %>%
  layer_dense(units = 4, activation = "tanh") %>%
  layer_dense(units = 2, activation = "tanh") %>%
  layer_dense(units = 4, activation = "tanh") %>%
  layer_dense(units = 8, activation = "tanh") %>%
  layer_dense(units = 15, activation = "tanh") %>%
  layer_dense(units = 29, activation = "tanh") %>%
  layer_dense(units = 58, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

# compile model
encoder2 %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam"
)

# fit model
encoder2 %>% fit(
  x = x_train, 
  y = x_train, 
  epochs = 10000,
  verbose = 0
)

encoder7 <- keras_model_sequential()
encoder7 %>%
  layer_dense(units = 58, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 29, activation = "tanh") %>%
  layer_dense(units = 15, activation = "tanh") %>%
  layer_dense(units = 7, activation = "tanh") %>%
  layer_dense(units = 15, activation = "tanh") %>%
  layer_dense(units = 29, activation = "tanh") %>%
  layer_dense(units = 58, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

# compile model
encoder7 %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam"
)

# fit model
encoder7 %>% fit(
  x = x_train, 
  y = x_train, 
  epochs = 10000,
  verbose = 0
)

model.eval(encoder7,x_train,x_train,4,'encoder7')
encoder7.mse <- evaluate(encoder7, x_train, x_train)

encoder2out <- keras_model(inputs = encoder2$input, outputs = get_layer(encoder2, index=6)$output)
encoder2bottleneck <- predict(encoder2out, x_train)

pdf('encoder2cond.pdf')
ggplot(data.frame(bottleneck1 = encoder2bottleneck[,1], bottleneck2 = encoder2bottleneck[,2]), aes(x = bottleneck1, y = bottleneck2, col = groups$Condition)) + geom_point()
dev.off()
pdf('bottleneckpheno.pdf')
ggplot(data.frame(bottleneck1 = intermediate_output[,1], bottleneck2 = intermediate_output[,2]), aes(x = bottleneck1, y = bottleneck2, col = groups$Phenotype)) + geom_point()
dev.off()

dists <- as.matrix(dist(intermediate_output))
neighbors <- apply(dists,2,order)

adj <- sapply(1:ncol(neighbors),function(i){
		      r <- dists[,i]
		      sel <- neighbors[-1:-7,i]
		      r[sel] <- 0
		      return(r)
})


g <- graph_from_adjacency_matrix(adj,'directed',T,F)
knn(g)

source('clustplots.R')

enrichCond(groups$Condition,adj,'bottleneck')


# set model
classifier <- keras_model_sequential()
classifier %>%
  layer_dense(units = 500, activation = "relu", input_shape = ncol(x_train)) %>%
  layer_dense(units = 100, activation = "relu") %>%
  layer_dense(units = 50, activation = "relu") %>%
  layer_dense(units = length(unique(clusts$leiden1.0)), activation = "softmax")

# compile model
classifier %>% compile(
  loss = "categorical_crossentropy", 
  optimizer = "adam",metrics='accuracy'
)

trainsel <- sample(1:nrow(x_train),nrow(x_train) %/% 10 * 9)
train <- x_train[trainsel,]
test <- x_train[-trainsel,]

clusts <- read.csv('2022_03_10/raw/dat/obs.csv'row.names=1)
clustenc <- sapply(0:max(clusts$leiden1.0), function(x) clusts$leiden1.0==x)
trainclust <- clusts$leiden1.0[trainsel]
testclust <- clusts$leiden1.0[-trainsel]

# fit model
classifier %>% fit(
  x = train, 
  y = trainclust, 
  epochs = 2000,
  verbose = 0
)

err <- read.delim('2022-09-07/encode2/log2error.txt')
nclust <- sapply(clusts2[2:11],sapply,max)

dir.pdf('log2error')
Heatmap(err,cluster_rows = F,cluster_columns = F,name='-log2(error)')
dev.off()

dir.pdf('nclust')
Heatmap(nclust,cluster_rows = F,cluster_columns = F,name='number of clusters')
dev.off()

dir.pdf('err_nclust')
Heatmap(-log2((2^(-err))/nclust),cluster_rows = F,cluster_columns = F,name='log2(error*numer of clusters)')
dev.off()
