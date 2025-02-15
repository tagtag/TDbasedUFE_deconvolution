#suppose the follwing files exist in the current directory
# Visium data sets
# GSM8341551_sample1_matrix.mtx.gz  GSM8341553_sample3_matrix.mtx.gz
# GSM8341552_sample2_matrix.mtx.gz  GSM8341554_sample4_matrix.mtx.gz
# GSM8341551_sample1_features.tsv.gz  GSM8341553_sample3_features.tsv.gz
# GSM8341552_sample2_features.tsv.gz  GSM8341554_sample4_features.tsv.gz
# GSM8341551_sample1_barcodes.tsv.gz  GSM8341553_sample3_barcodes.tsv.gz
# GSM8341552_sample2_barcodes.tsv.gz  GSM8341554_sample4_barcodes.tsv.gz
# scRNA-seq
# GSM8341613_1_3DE_genematrix.csv.gz GSM8341614_2_3DE_genematrix.csv.gz
# GSM8341772_control_3DE_genematrix.csv.gz  GSM8341774_injected_3DE_genematrix.csv.gz

#Storing Visium data sets
files <- list.files("./",pattern="sample[1-4]_matrix")
X<-NULL
for (i in seq_along(files))
{
  cat(i," ")
  x <- read.csv(files[i],sep="",comment.char="%",header=F)
  xa <- read.csv( gsub("_matrix.mtx.gz","_barcodes.tsv.gz",files[i]),header=F)
  xb <- read.csv(gsub("_matrix.mtx.gz","_features.tsv.gz",files[i]),sep="\t",header=F)
  y <-  read.csv( gsub("_matrix.mtx.gz","_tissue_positions_list.csv.gz",files[i]),header=F)
  x0<- spMatrix(nrow=dim(xb)[1],ncol=dim(y)[1],i=x[,1],j=match(xa[x[,2],],y[,1]),x=x[,3])
  x0 <- x0/mean(x0)
  X <- cbind(X,x0)
}
save(file="X",X)

#storing  scRN-seq
xx <- read.csv("GSM8341613_1_3DE_genematrix.csv.gz",sep=",");save(file="xx1",xx)
xx <- read.csv("GSM8341614_2_3DE_genematrix.csv.gz",sep=",");save(file="xx2",xx)
xx <- read.csv("GSM8341772_control_3DE_genematrix.csv.gz",sep=",");save(file="xx3",xx)
xx <- read.csv("GSM8341774_injected_3DE_genematrix.csv.gz",sep=",");save(file="xx4",xx)


require(irlba)
SVD <- irlba(X)

#plot v_{2j''} 
pdf(file="SVD.pdf",width=7.5,height=5)
m <- matrix(c(
  1, 2, 5,
  3, 4, 6
), nrow = 2, ncol = 3, byrow = TRUE)

# layout set up
layout(m)

for (i in seq_along(files))
{
  y <-  read.csv( gsub("_matrix.mtx.gz","_tissue_positions_list.csv.gz",files[i]),header=F)
  a <- SVD$v[4992*(i-1)+c(1:4992),2]
  CUT <- cut(a,10)
  b<- cbind(levels(CUT),rainbow(length(levels(CUT))))
  col <- b[match(CUT,b[,1]),2]
  plot(y[,6],rev(y[,5]),col=col,pch=16,cex=0.5,xlab="x",ylab="y",cex.lab=1.5)
}
plot(rep(0,10),1:10,col=b[,2],pch=16,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,0.5),bty="n",ylim=c(0.5,10.5))
text(rep(0,10),1:10,b[,1],pos=4)
dev.off()

#selecting 277 genes
k0<-2
th <- function(sd){
  P2<- pchisq(((SVD$u[,k0]-mean(SVD$u[,k0]))/sd)^2,1,lower.tail=F)
  hc<- hist(1-P2,breaks=100,plot=F)
  return(sd(hc$count[1:sum(hc$mid<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.001,th)$par
P1<- pchisq(((SVD$u[,k0]-mean(SVD$u[,k0]))/sd)^2,1,lower.tail=F)
table(p.adjust(P1,"BH")<0.01)
#FALSE  TRUE
#30776   277


#plot v^{sel}_{1j''} and  v^{sel}_{4j''} as well as Spearman correlation coefficients 
#choose id<-1 or id <-4 
SVD1 <- irlba(X[p.adjust(P1,"BH")<0.01,])
#id <-1
id<-4
a <- SVD1$v[,id]
CUT <- cut(a,10)
b<- cbind(levels(CUT),rainbow(length(levels(CUT))))
col <- b[match(CUT,b[,1]),2]
#pdf(file="SVD_sel.pdf",width=7.5,height=5) #id <-1
pdf(file="SVD_sel_4.pdf",width=7.5,height=5) #id <-4
#par(mfrow=c(2,2))
m <- matrix(c(
  #  1, 2, 3, #id<-1
  #  4, 5, 6  #id<-1
  1, 2, 5, #id<-4
  3, 4, 6  #id<-4
), nrow = 2, ncol = 3, byrow = TRUE)

# layout set up
layout(m)
files1 <- list.files("./",pattern="GSM8341[6-7]")
for (i in seq_along(files1))
{
  cat("start reading file")
  #xx <- read.csv(files1[i],sep=",")
  #save(file="xx1",xx)
  load(paste("xx",i,sep=""))
  cat("end reading file")
  xb <- read.csv(gsub("_matrix.mtx.gz","_features.tsv.gz",files[i]),sep="\t",header=F)
  xx <- xx[match(xb[p.adjust(P1,"BH")<0.01,1],xx[,1]),]
  xx1 <- xx[,-1]
  xx1 <- data.matrix(xx1)
  X1 <- X[p.adjust(P1,"BH")<0.01,4992*(i-1)+c(1:4992)]
  index <- colSums(X1)!=0
  
  svdY   <- irlba(xx1,5)       
  Uprime <- svdY$u  
  y <-  read.csv( gsub("_matrix.mtx.gz","_tissue_positions_list.csv.gz",files[i]),header=F)
  col <- b[match(CUT[4992*(i-1)+c(1:4992)],b[,1]),2]
  plot(y[,6],rev(y[,5]),col=col,pch=16,cex=0.5,xlab="x",ylab="y",cex.lab=1.5)
  LM <- lm(data.matrix(X1[,index])~Uprime)
  SLM <- summary(LM)
  corP<-NULL
  for (j in seq_along(SLM))
  {
    #cat(j," ")
    COR <- cor.test(rank(LM$fitted[,j]),rank(X1[,index][,j])) 
    corP <- cbind(corP,c(COR$estimate,COR$p.value))
  }
  
  #hist(abs(corP[1,]),breaks=100,xlab="|cor|",cex.lab=1.5) #id <-1
  #plot(rep(0,10),1:10,col=b[,2],pch=16,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,0.5),bty="n",ylim=c(0.5,10.5)) #id <-1
  #text(rep(0,10),1:10,b[,1],pos=4) #id<-1
  print(table(p.adjust(corP[2,],"BH")<0.05))
}
plot(rep(0,10),1:10,col=b[,2],pch=16,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(0,0.5),bty="n",ylim=c(0.5,10.5)) #id<-4
text(rep(0,10),1:10,b[,1],pos=4) #id<-4
dev.off()
#The numbef or spots associated with significant adjusted P-values
table(p.adjust(corP[2,],"BH")<0.05)
#k=1
#FALSE  TRUE
#971  2076
#k=2
#FALSE  TRUE
#1271  1769
#k=3
#FALSE  TRUE
#1261  1349
#k=4
#FALSE  TRUE
#1230  1351

#assigning cell types to individual cells in  scRN-seq
for (i in c(1:4))
{
load(paste("xx",i,sep=""))
library(celldex)
ref <-celldex::MouseRNAseqData()
require(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
ensembl_genes <- xx[,1]  # 例: Ensembl Gene Name
results <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),  
  filters = "ensembl_gene_id",
  values = ensembl_genes,
  mart = ensembl
)
genes <- results[match(xx[,1],results[,1]),2]
xx<- xx[!duplicated(genes),]
genes <- results[match(xx[,1],results[,1]),2]
xx<- xx[!is.na(genes),]
rownames(xx) <- results[match(xx[,1],results[,1]),2]
results <- SingleR(test = xx[,-1], ref = ref, labels = ref$label.main)
save(file=paste("results",i,sep=""),results)]
}

#comuting and plotting cell fraction
cell_frac<- rep(list(NA),4)
for (i in seq_along(files1))
{
  cat("start reading file")
  load(paste("xx",i,sep=""))
  cat("end reading file")
  xb <- read.csv(gsub("_matrix.mtx.gz","_features.tsv.gz",files[i]),sep="\t",header=F)
  xx <- xx[match(xb[p.adjust(P1,"BH")<0.01,1],xx[,1]),]
  xx1 <- xx[,-1]
  xx1 <- data.matrix(xx1)
  X1 <- X[p.adjust(P1,"BH")<0.01,4992*(i-1)+c(1:4992)]
  index <- colSums(X1)!=0
  
  svdY   <- irlba(xx1,5)       
  Uprime <- svdY$u  
  LM <- lm(data.matrix(X1[,index])~Uprime)
  SLM <- summary(LM)
  corP<-NULL
  for (j in seq_along(SLM))
  {
    COR <- cor.test(rank(LM$fitted[,j]),rank(X1[,index][,j])) 
    corP <- cbind(corP,c(COR$estimate,COR$p.value))
  }
 g_cont <- svdY$v%*%LM$coefficients[-1,]+matrix(LM$coefficients[1,],nrow= dim(svdY$v)[1],ncol=dim(LM$coefficients[-1,])[2],byrow=T)
  g_cont<- g_cont[,p.adjust(corP[2,],"BH")<0.05]
  load(paste("results",i,sep=""))
  cell_frac[[i]] <- apply(g_cont,2,function(x){unlist(lapply(split(x,results$pruned.labels),sum))})
}
save(file="cell_frac",cell_frac)
N1 <-100
i<-1;Y <- t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean))
i<-2; Y <- merge(Y,t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean)),by="row.names",all=T)
rownames(Y) <- Y[,1]
Y <- Y[,-1]
i<-3; Y <- merge(Y,t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean)),by="row.names",all=T)
rownames(Y) <- Y[,1]
Y <- Y[,-1]
i<-4; Y <- merge(Y,t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean)),by="row.names",all=T)
rownames(Y) <- Y[,1]
Y <- Y[,-1]
coln <- length(rownames(Y))
cell_names <- rownames(Y)
Y <- apply(Y,2,as.numeric)
Y[is.na(Y)]<-0
Ymax <- max(colSums(Y),na.rm=T)*3
colnames(Y)<-rep(1:100,4)
pdf(file="cell_frac.pdf")
for(i in c(1:length(cell_frac)))
{
  barplot(Y[,(i-1)*100+c(1:100)],col=rainbow(coln),ylim=c(0,Ymax))
  legend(20,Ymax*0.9,cell_names,pch=16,col=rainbow(coln))
}
dev.off()

#RCTD 
require(spacexr)
require(Matrix)
require(doParallel)
load("X")
i<-1
#i<-2
#i<-3
#i<-4
counts <- X[,(i-1)*4992+c(1:4992)]
coords <-  read.csv( gsub("_matrix.mtx.gz","_tissue_positions_list.csv.gz",files[i]),header=F)
rownames(coords) <- coords[,1]; coords[,1] <-NULL
coords <- coords[,-c(1:3)]
xb <- read.csv(gsub("_matrix.mtx.gz","_features.tsv.gz",files[i]),sep="\t",header=F)
xa <- read.csv( gsub("_matrix.mtx.gz","_barcodes.tsv.gz",files[i]),header=F)
rownames(counts) <- xb[,1]
index2 <- colSums(counts)!=0
counts <- counts[,index2]
coords <- coords[index2,]
rownames(counts) <- xb[,1]
colnames(counts) <- xa[,1]
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
puck <- SpatialRNA(coords, counts, nUMI,require_int=F)
load(paste("xx",i,sep=""))
load(paste("results",i,sep=""))
rownames(xx) <- xx[,1]; xx[,1] <- NULL
cell_types <-data.frame(results@rownames,results$labels)
cell_types <- setNames(cell_types[[2]], cell_types[[1]])
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- data.frame(colnames(xx),colSums(xx))
nUMI <- setNames(nUMI[[2]], nUMI[[1]])
reference <- Reference(xx, cell_types, nUMI)
#myRCTD <- create.RCTD(puck, reference, max_cores = 2, CELL_MIN_INSTANCE=1,cell_type_names=levels(reference@cell_types)[-3]) #i<-1
#myRCTD <- create.RCTD(puck, reference, max_cores = 2, CELL_MIN_INSTANCE=1,cell_type_names=levels(reference@cell_types)[-c(7:9)]) #i<-2
#myRCTD <- create.RCTD(puck, reference, max_cores = 2, CELL_MIN_INSTANCE=1,cell_type_names=levels(reference@cell_types)[-c(2,3,4,7,11)]) #i<-3
myRCTD <- create.RCTD(puck, reference, max_cores = 2, CELL_MIN_INSTANCE=1,cell_type_names=levels(reference@cell_types)[-c(2,3)]) #i<-4
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
save(file=paste("myRCTD",i,sep=""),myRCTD)
i<-1;Y <- t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean))
i<-2; Y <- merge(Y,t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean)),by="row.names",all=T)
rownames(Y) <- Y[,1]
Y <- Y[,-1]
i<-3; Y <- merge(Y,t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean)),by="row.names",all=T)
rownames(Y) <- Y[,1]
Y <- Y[,-1]
i<-4; Y <- merge(Y,t(t(cell_frac[[i]][,1:N1])*apply(apply(cell_frac[[i]][,1:N1],2,sign),2,mean)),by="row.names",all=T)
rownames(Y) <- Y[,1]
Y <- Y[,-1]
coln <- length(rownames(Y))
cell_names <- rownames(Y)
Y <- apply(Y,2,as.numeric)
Y[is.na(Y)]<-0
Ymax <- max(colSums(Y),na.rm=T)*3
colnames(Y)<-rep(1:100,4)
pdf(file="cell_frac.pdf")
for(i in c(1:length(cell_frac)))
{
  barplot(Y[,(i-1)*100+c(1:100)],col=rainbow(coln),ylim=c(0,Ymax))
  legend(20,Ymax*0.9,cell_names,pch=16,col=rainbow(coln))
}
dev.off()

#SPOTlight
require(scran)
require(SPOTlight)
for (i in c(1:4))
{
load(paste("xx",i,sep=""))
load(paste("results",i,sep=""))
rownames(xx) <- xx[,1]; xx[,1] <- NULL
mgs <- scoreMarkers(xx,results$labels)
load("X")
counts <- X[,(i-1)*4992+c(1:4992)]
xb <- read.csv(gsub("_matrix.mtx.gz","_features.tsv.gz",files[i]),sep="\t",header=F)
xa <- read.csv( gsub("_matrix.mtx.gz","_barcodes.tsv.gz",files[i]),header=F)
rownames(counts) <- xb[,1]
index2 <- colSums(counts)!=0
counts <- counts[,index2]
rownames(counts) <- xb[,1]
colnames(counts) <- xa[,1]
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
res <- SPOTlight(
  x = data.matrix(xx),
  y = data.matrix(counts),
  groups = results$labels,
  mgs = mgs_df,
  hvg = NULL,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")
save(file=paste("res",i,sep=""),res)
}

i<-1
load(paste("res",i,sep=""))
Y <-  data.frame(t(res$mat))[,1:100]
for (i in c(2:4))
{
  load(paste("res",i,sep=""))
  Y <- merge(Y,data.frame(t(res$mat))[,1:100],by="row.names",all=T)
  rownames(Y) <- Y[,1]
  Y <- Y[,-1]
}
coln <- length(rownames(Y))
cell_names <- rownames(Y)
Y <- apply(Y,2,as.numeric)
Y[is.na(Y)]<-0
Ymax <- max(colSums(Y),na.rm=T)*3
colnames(Y)<-rep(1:100,4)

pdf(file="cell_frac_SPOTlight.pdf")
for(i in c(1:length(files)))
{
  barplot(data.matrix(Y[,(i-1)*100+c(1:100)]),col=rainbow(coln),ylim=c(0,Ymax),xaxt="n")
  legend(20,Ymax,cell_names,pch=16,col=rainbow(coln))
}
dev.off()

#SpaCET
#i<-1
#i<-2
#i<-3
i<-4
require(SpaCET)
load(paste("xx",i,sep=""))
rownames(xx)<- xx[,1]
xx <- xx[,-1]
xx <- data.matrix(xx)
load(paste("results",i,sep=""))
sc_annotation <- data.frame(rownames(results),results[,2])
celltype <- names(table(sc_annotation[,2]))
#i=1
#sc_lineageTree <- list(Astrocytes="Astrocytes","B cells"= "B cells","Epithelial cells"="Epithelial cells",
#                       "Fibroblasts"="Fibroblasts","Macrophages"="Macrophages","Microglia"="Microglia","Neurons"="Neurons",
#                       "NK cells"="NK cells","Oligodendrocytes"="Oligodendrocytes","T cells"="T cells")
#i=2
#sc_lineageTree <- list(Astrocytes="Astrocytes","B cells"= "B cells","Monocytes"="Monocytes",
#                       "Fibroblasts"="Fibroblasts","Macrophages"="Macrophages","Microglia"="Microglia","Neurons"="Neurons",
#                       "NK cells"="NK cells","Oligodendrocytes"="Oligodendrocytes","T cells"="T cells","Granulocytes"="Granulocytes")
#i=3
#sc_lineageTree <- list(Astrocytes="Astrocytes","Monocytes"="Monocytes",
#                       "Fibroblasts"="Fibroblasts","Macrophages"="Macrophages","Microglia"="Microglia","Neurons"="Neurons",
#                       "NK cells"="NK cells","Oligodendrocytes"="Oligodendrocytes","T cells"="T cells","Granulocytes"="Granulocytes",
#                       "Dendritic cells"="Dendritic cells", "Epithelial cells"="Epithelial cells","Endothelial cells"="Endothelial cells")
#i=4
sc_lineageTree <- list(Astrocytes="Astrocytes", "B cells"= "B cells",
                       "Fibroblasts"="Fibroblasts","Microglia"="Microglia","Neurons"="Neurons",
                       "Oligodendrocytes"="Oligodendrocytes","T cells"="T cells")

load("X")
X <- X[,4992*(i-1)+c(1:4992)]
feature <- read.csv(gsub("_matrix.mtx.gz","_features.tsv.gz",files[i]),sep="\t",header=F)
rownames(X) <- feature[,1]
X <- X[,colSums(X)>0]
barcode <-  read.csv( gsub("_matrix.mtx.gz","_barcodes.tsv.gz",files[i]),header=F)
colnames(X) <- barcode[,1]
location<- read.csv( gsub("_matrix.mtx.gz","_tissue_positions_list.csv.gz",files[i]),header=F)
rownames(location) <- location[,1]
location <- location[match(colnames(X),location[,1]),5:6]

SpaCET_obj <- create.SpaCET.object(
  counts=X,
  spotCoordinates=location,
  imagePath=NA,
  platform = "oldST"
)

colnames(sc_annotation) <- c("cellID","bio_celltype")
SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
  SpaCET_obj, 
  sc_counts=xx, 
  sc_annotation=sc_annotation, 
  sc_lineageTree=sc_lineageTree, 
  coreNo=8
)

#save(file="SpaCET_obj1",SpaCET_obj)
#save(file="SpaCET_obj2",SpaCET_obj)
#save(file="SpaCET_obj3",SpaCET_obj)
save(file="SpaCET_obj4",SpaCET_obj)

i<-1
load(paste("SpaCET_obj",i,sep=""))
Y <-  SpaCET_obj@results[[1]]$propMat[,1:100]
for (i in c(2:4))
{
  load(paste("SpaCET_obj",i,sep=""))
  Y <- merge(Y,data.frame(SpaCET_obj@results[[1]]$propMat)[,1:100],by="row.names",all=T)
  rownames(Y) <- Y[,1]
  Y <- Y[,-1]
}
coln <- length(rownames(Y))
cell_names <- rownames(Y)
Y <- apply(Y,2,as.numeric)
Y[is.na(Y)]<-0
Ymax <- max(colSums(Y),na.rm=T)*3
colnames(Y)<-rep(1:100,4)

pdf(file="cell_frac_SpaCET.pdf")
for(i in c(1:length(files)))
{
  barplot(data.matrix(Y[,(i-1)*100+c(1:100)]),col=rainbow(coln),ylim=c(0,Ymax),xaxt="n")
  legend(20,Ymax,cell_names,pch=16,col=rainbow(coln))
}
dev.off()

#cell2location

i<-1
x <- read.csv(paste("cell_abundances",i,".csv",sep=""))
Y <- data.matrix(t(x[2:101,]))
Y <- Y[-1,]
for (i in c(2:4))
{
  x <- read.csv(paste("cell_abundances",i,".csv",sep=""))
  Y <- 
    Y <- Y[-1,]
  load(paste("SpaCET_obj",i,sep=""))
  Y <- merge(Y,data.matrix(t(x[2:101,]))[-1,][,1:100],by="row.names",all=T)
  rownames(Y) <- Y[,1]
  Y <- Y[,-1]
}
coln <- length(rownames(Y))
cell_names <- rownames(Y)
Y <- apply(Y,2,as.numeric)
Y[is.na(Y)]<-0
Ymax <- max(colSums(Y),na.rm=T)*3
colnames(Y)<-rep(1:100,4)

pdf(file="cell_frac_cell2location.pdf")
for(i in c(1:length(files)))
{
  celltypes<- as.vector(data.frame(strsplit(cell_names,"_"))[10,])
  barplot(data.matrix(Y[,(i-1)*100+c(1:100)]),col=rainbow(coln),ylim=c(0,Ymax),xaxt="n")
  legend(20,Ymax,rev(celltypes) ,pch=16,col=rev(rainbow(coln)))
}
dev.off()

