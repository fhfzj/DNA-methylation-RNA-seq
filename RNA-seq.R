
 ############################################
 #                                          #
 #  DEFFERENTIAL EXPRESION ANALYSIS CHAGAS  #            
 #                                          #
 ############################################



setwd("/mnt/beegfs/fheredia/datos_expresion/")

library(limma)
library(edgeR) # version 4.0.1
library(gtools)
library(tweeDEseqCountData)
library(RColorBrewer)
library(enrichR)

 
################################
### CARGAR TARGETS Y READCOUNTS
###############################
 
targets <- read.delim("targets_sel_nombre_met.txt", sep = "\t")
targets

#readcount generado en miARma
readcount<-read.delim("/mnt/beegfs/fheredia/datos_expresion/miARma_results/Readcount_results/str-ReadCount.tab")

#MODIFICAR DATOS DE LOS TARGETS Y READCOUNTS
#Ordenar columnas en readcounts
#mixedosrt: ordena nombres con letras y numeros -> library(gtools)
readcount_ordenado <- readcount[, mixedsort(colnames(readcount))]

#Modificar nombres en targets - para que coincida con los nombres del readcount 
# Función para modificar los nombres
modificar_nombres <- function(nombre) {
  # Eliminar las dos últimas letras y añadir "nat"
  nuevo_nombre <- paste0(substr(nombre, 1, nchar(nombre) - 2), "nat")
  return(nuevo_nombre)
}

# Crear una nueva columna con los nombres modificados
targets$Name_2 <- sapply(targets$Name, modificar_nombres)

#hacemos nueva columna con CCC y healthy
targets$Status <- ifelse(targets$Type %in% "HD", targets$Type, "CCC")

#readcount con los indiviudos que estan en targets
readcount_def<-readcount_ordenado[,which(colnames(readcount_ordenado)%in%targets$Name_2)]


#########################
##ANOTACION y DGE OBJECT 
##########################

#Cargar ANOTACION de miARma
size<-read.delim("/mnt/beegfs/fheredia/datos_expresion/miARma_results/Readcount_results/str-Size.tab")
m<-match(rownames(readcount_def),size$Gene)
annot<-size[m,]
annot[is.na(annot$Length),"Length"]<-0

#CREAR DGE OBJETC
#Obtener CPMs
counts_CPM <- cpm(readcount_def)
y <- DGEList(counts = counts_CPM, genes = annot, group = targets$Type)


#################################
### FILTRERING AND NORMALIZATION 
#################################

#Filtro para eliminar genes con menos de 1 cpm en el grupo con menor muestras
keep <- filterByExpr(y, min.count=1) #para mas informacion sobre la funcion mira ?filterByExpr y https://f1000research.com/articles/5-1438
# "The function accesses the group factor contained in y in order to compute the minimum group size, but the filtering is performed independently of which sample belongs to which group so that no bias is introduced."
table(keep) #nos quedamos con 20870 genes
y_raw <- y[keep, , keep.lib.sizes=FALSE] #Recomputar el tamaño de la librería

#NORMALIZACION 
y <- calcNormFactors(y_raw)


####################
### QUALITY CONTROL    
####################

#PASAR GRÁFICOS A PDF
setwd("/mnt/beegfs/fheredia/DEAnalysis/")
pdf(file=paste("graficosQC_DEA", "_", format(Sys.Date(), "%d.%m.%Y"), ".pdf", sep = ""), width=7, height=6, paper="a4")

#TAMAÑO DE LA LIBRERÍA
#si el nº de reads de cada muestra es muy diferente entre sí hay que hacer downsampling
lib <- y$samples$lib.size
barplot(lib, names.arg = targets$Name, main = "Tamaño de la librería", ylab = "Reads", xlab = "", col = "steelblue", las = 2, cex.names = 0.8, cex.axis = 0.6, cex.main=0.9)

#BOXPLOT
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
#Boxplot of the unnormalised samples (cpm)
log_CPM_raw<-cpm(y_raw, log=TRUE, )
boxplot(log_CPM_raw, col=pal[factor(targets$Status)], main="Boxplot of unnormalised samples",las=2, names=targets$Name, ylab="logCPM", xlab="", cex.axis=0.6)
#Boxplot of the normalized samples
logCPM<-cpm(y, log=TRUE)
boxplot(logCPM, col=pal[factor(targets$Status)], main="Boxplot of normalized samples", las=2, names=targets$Name, ylab="logCPM", xlab="", cex.axis=0.6)

par(mfrow=c(1,1))

#DENDROGRAMA
clust<- hclust(na.omit(dist(t(logCPM))))
plot(clust, main="HCLUST logCPM", labels=targets$Name)

#PCA
data.PC = prcomp(t(as.matrix(logCPM)))
plot.default(data.PC$x,col=pal[factor(targets$Type)],main="PCA logCPM", pch=16)
text(data.PC$x[,1],data.PC$x[,2],labels=targets$Type, cex=0.7, pos=1)

#PCOA
plotMDS(logCPM, gene.selection="common", 
        col=pal[factor(targets$Type)], cex=0.8, main="PCOA logCPM")
legend("top", legend=levels(factor(targets$Type)), text.col=pal,
       cex=0.65, bg="white")

#cerrar pdf
dev.off()


#######################
###LINEAR MODELING
#######################

#K1,K2,K3 vs HD
#design matrix
design <- model.matrix(~0+Type+Sexo+Edad, data=targets)
colnames(design) <- c("H", "K1", "K2", "K3", "M","Age")
design

#"Use voom to convert the read counts to log2-cpm, with associated weights, ready for linear modelling"
# voom tmb se puede usar con los counts directamente o con los cpm [viñeta limma P.72]
v <- voom(y, design, plot = TRUE)
#¿se necesitaría un estudio de correlación?
fit <- lmFit(v, design)

contMatrix <- makeContrasts(K1_vs_H=K1-H,
                            K2_vs_H=K2-H,
                            K3_vs_H=K3-H,
                            levels=design)
contMatrix
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))
DEG_K1vsH<-topTable(fit2,coef="K1_vs_H",sort.by="p", num=Inf)
  DEG_K1vsH$FDR <- p.adjust(DEG_K1vsH$P.Value, method = 'fdr')
DEG_K2vsH<-topTable(fit2,coef="K2_vs_H",sort.by="p", num=Inf)
  DEG_K2vsH$FDR <- p.adjust(DEG_K2vsH$P.Value, method = 'fdr')
DEG_K3vsH<-topTable(fit2,coef="K3_vs_H",sort.by="p", num=Inf)
  DEG_K3vsH$FDR <- p.adjust(DEG_K3vsH$P.Value, method = 'fdr')

###CCC vs H 
design_2 <- model.matrix(~0+Status+Sexo+Edad, data=targets)
colnames(design_2) <- c("CCC", "H", "M","Age")
design_2
v_2 <- voom(y, design_2, plot = TRUE)
fit_2 <- lmFit(v_2, design_2)
contMatrix_2 <- makeContrasts(CCC_vs_H=CCC-H,
                            levels=design_2)
contMatrix_2
fit2_2 <- contrasts.fit(fit_2, contMatrix_2)
fit2_2 <- eBayes(fit2_2)
summary(decideTests(fit2_2))
DEG_CCCvsH<-topTable(fit2_2,coef="CCC_vs_H",sort.by="p", num=Inf)
DEG_CCCvsH$FDR <- p.adjust(DEG_CCCvsH$P.Value, method = 'fdr')

###K2K3 vs H
design_3 <- model.matrix(~0+Join+Sexo+Edad, data=targets)
colnames(design_3) <- c("H", "K1", "K2K3","M", "Age")
design_3
v_3 <- voom(y, design_3, plot = TRUE)
fit_3 <- lmFit(v_3, design_3)
contMatrix_3 <- makeContrasts(K2K3_vs_H=K2K3-H,
                              levels=design_3)
contMatrix_3
fit2_3 <- contrasts.fit(fit_3, contMatrix_3)
fit2_3 <- eBayes(fit2_3)
summary(decideTests(fit2_3))
DEG_K2K3vsH<-topTable(fit2_3,coef="K2K3_vs_H",sort.by="p", num=Inf)
DEG_K2K3vsH$FDR <- p.adjust(DEG_K2K3vsH$P.Value, method = 'fdr')


###################
###ENRIQUECIMIENTO
###################


#to be continued...
