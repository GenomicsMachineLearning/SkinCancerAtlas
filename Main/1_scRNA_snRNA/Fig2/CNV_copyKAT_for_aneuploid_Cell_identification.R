setwd("/QRISdata/Q4386/skin_atlas/SCC_BCC/copykat/")
library(copykat)
g<-read.csv("gene_mat.txt", sep="\t", header=TRUE)
rownames(g)<-g$GENE
g<-g[,-1]
copykat.test <- copykat(rawmat=g, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=1) #hg20 built-in copkat is the hg38 coords
