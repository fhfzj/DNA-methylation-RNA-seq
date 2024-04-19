
 ######################################
 #                                    #
 #  ANÁLISIS METILACIÓN ADN CHAGAS    #
 #                                    #
 ######################################


setwd("/mnt/beegfs/fheredia")
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(openxlsx)
library(ggplot2)
library(ggfortify)


################
###CARGAR DATOS 
################
 
#Hay datos de 2019 y 2021 en carpetas separadas
dataDirectory_2019<-"/mnt/beegfs/fheredia/datos_metilacion/2019-123-ILL_idats"
#list.files(dataDirectory_2019, recursive=TRUE)
targets_2019 <- read.metharray.sheet(dataDirectory_2019, pattern="Targets2.csv")
dataDirectory_2021<-"/mnt/beegfs/fheredia/datos_metilacion/2021-176-ILL_METGRA_N=116"
#list.files(dataDirectory_2021, recursive=TRUE)
targets_2021 <- read.metharray.sheet(dataDirectory_2021, pattern="Targets3.csv")
#unificar los 2 targets en uno:  
targets<-rbind(targets_2019, targets_2021)

 
rgSet <- read.metharray.exp(targets=targets, force=TRUE) 
#calculo del valor p
detP<-detectionP(rgSet)


#############################
### FILTRADO Y NORMALIZACIÓN
#############################

#eliminar muestras de mala calidad:
keep <- colMeans(detP) < 0.05   
table(keep)   #no se elimina ninguna, con p<0.01 tampoco
rgSet<-rgSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]

#NORMALISATION
#Doble normalizacion Noob-Quantile
rgSet_Noob<-preprocessNoob(rgSet)
mSetSq<-preprocessQuantile(rgSet_Noob, fixOutliers = TRUE,
                   removeBadSamples = TRUE,
                   badSampleCutoff = 10.5,
                   quantileNormalize = TRUE,
                   stratified = TRUE,
                   mergeManifest = FALSE,
                   sex = NULL)
#"if gender is unspecified (NULL), a guess is made using by the getSex function using copy number information"

#PREDICCION DE SEXO
predictedSex <- getSex(mSetSq, cutoff = -2)$predictedSex
table(targets$Gender==predictedSex) #NO COINCIDE EN UN INDIVIDUO 
#lo introduzco en targets
targets$predictedSex<-predictedSex
#este será el sexo que se use en las design matrix

#FILTRADO DE SONDAS
#poner mismo orden las sondas en el mSetSq y en detP
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
#sondas que no han fallado en ningun caso
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep) #7920 sondas eliminadas
mSetSqFlt <- mSetSq[keep,]

#eliminar sondas que hibridan en los cromosomas sexuales 
#informacion sobre las sondas usadas
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])
table(keep) #se eliminan 10183 sondas
mSetSqFlt <- mSetSqFlt[keep,]

#eliminar sondas con SNPs
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

#eliminar sondas cross reactive 
xReactiveProbes <- read.csv(file=paste(dataDirectory_2019,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)   #se eliminan 24660 sondas
mSetSqFlt <- mSetSqFlt[keep,] 

#Calculo valor M y beta de los datos filtrados
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)


#####################
### QUALITY CONTROL
#####################

#qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$Status, 
#         pdf="qcReportFH1.pdf")

#PASAR GRÁFICOS A PDF
setwd("/mnt/beegfs/fheredia/Scripts_metilacion/")
pdf(file=paste("graficosQC_metilacion", "_", format(Sys.Date(), "%d.%m.%Y"), ".pdf", sep = ""), width=10, height=6, paper="a4")

#ESTUDIO DE LA NORMALIZACIÓN
#para comparar el antes y despues de la normalizacion obtenemos un metylset de la raw data:
mSetRaw <- preprocessRaw(rgSet)
#Representación la densidad del valor beta antes y después
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Status,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Status)), 
       text.col=brewer.pal(8,"Dark2"),cex = 0.5)
densityPlot(getBeta(mSetSq), sampGroups=targets$Status,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Status)), 
       text.col=brewer.pal(8,"Dark2"), cex = 0.5)
par(mfrow=c(1,1))

#DENSIDAD DE LOS VALORES M Y BETA 
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Status, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Status)), 
       text.col=brewer.pal(8,"Dark2"), cex = 0.5)
densityPlot(mVals, sampGroups=targets$Status, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Status)), 
       text.col=brewer.pal(8,"Dark2"), cex = 0.5)
par(mfrow=c(1,1))

#PCA
data.PC = prcomp(t(as.matrix(mVals)))
plot.default(data.PC$x,col=pal[factor(targets$Status)],main="PCA-Status mVals", pch=16)
text(data.PC$x[,1],data.PC$x[,2],labels=targets$Status, cex=0.7, pos=1)

#PCA - otra representación
mvals_df <- as.data.frame(mVals)
colnames(mvals_df) <- targets$Type
data.t<-as.data.frame(t(na.omit(mvals_df)))
rownames(data.t)<-paste(targets$Sample_ID)
data.PC = prcomp(data.t)
data.t$Status<-targets$Status
autoplot(data.PC,label=T,data=data.t,colour='Status',label.size=3) + 
  ggtitle(label = "PCA - Status")

#PCOA
pal <- brewer.pal(8,"Dark2")
plotMDS(mVals, top=1000, gene.selection="common", 
        col=pal[factor(targets$Status)], cex=1, main="PCOA mVals")
legend("right", legend=levels(factor(targets$Status)), text.col=pal,
       cex=0.6, bg="white")

#DENDROGRAMA CON VALOR M (siempre hacer con valor M)
df_mVals <- as.data.frame(mVals)
clust<- hclust(na.omit(dist(t(df_mVals))))
plot(clust, main="HCLUST", labels=targets$Status)
#dendrograma con log2    
df_mVals_log2 <- log2(df_mVals)
clust_log2<- hclust(na.omit(dist(t(df_mVals_log2))))
plot(clust_log2, main="HCLUST_log2", labels=targets$Status)

#cerrar pdf
dev.off()


#########################
### SENTITIVITY TESTS
#########################

#estudio de las covariables que se deberían de tener en cuenta en el analisis
pdata<-targets
pval_var <- function(input = bVals, pdata = pdata, n_pca = 3, method = c("bonferroni", "holm", "hochberg", "BH", "BY", "fdr", "none")){
  
  pcoordinates <- mxec.dreduction(input)
  batch_effects <- mxec.batcheffects(pcoordinates, pdata)
  pval_final <- as.data.frame(batch_effects$pc$pca$pvalues)[,1:n_pca]
  colnames(pval_final) <- paste("pca", colnames(pval_final), sep="")
  FDR <- apply(pval_final, 2, function(x) p.adjust(x, method = method))
  matrix <- as.list(c())
  matrix$pval <- pval_final
  matrix$fdr <- FDR
  matrix$pval = as.data.frame(matrix$pval)
  matrix$fdr = as.data.frame(matrix$fdr)
  summary <- as.data.frame(summary(prcomp(na.omit(input)))$importance)
  pc1 <- summary[2,1]*100
  pc2 <- summary[2,2]*100
  pc3 <- summary[2,3]*100
  colnames(matrix$pval) <- c(paste("PC1", paste0(pc1,"%"), sep=" "), paste("PC2", paste0(pc2,"%"), sep=" "),paste("PC3", paste0(pc3,"%"), sep=" "))
  colnames(matrix$fdr) <- c(paste("PC1", paste0(pc1,"%"), sep=" "), paste("PC2", paste0(pc2,"%"), sep=" "),paste("PC3", paste0(pc3,"%"), sep=" "))
  matrix
}

# Heatmap tes de sensibilidad

var_heatmap <- function(input = matrix, title = "final list", pval = F){
  library(ggplot2)
  if(isTRUE(pval)){cormat <- input$pval}else{cormat <- input$fdr}
  cormat2 <- -log10(cormat)
  cormat$id =rownames(cormat)
  melted_cormat <- reshape2::melt(cormat, na.rm = TRUE)
  melted_cormat2 <- reshape2::melt(cormat2, na.rm = TRUE)
  melted_cormat$logfdr <- melted_cormat2[,2]
  if(isTRUE(pval)){colnames(melted_cormat) <- c('Var1','Var2','pval', "logpval")}else{
    colnames(melted_cormat) <- c('Var1','Var2','fdr', "logfdr")
  }
  
  
  buylrd <- c("#ffffff", "#f5c4c4", "#f59a9a", "#fa0000")
  colors.martin <- colorRampPalette(buylrd)(30)
  
  if(isTRUE(pval)){
    ggplot(melted_cormat, aes(Var2, Var1, fill = logpval)) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colours = colors.martin) +
      geom_text(aes(Var2, Var1, label = round(pval,3)), color = "black", size = 4) +
      xlab("") +
      ylab('') +
      ggtitle(title)
  }else{
    ggplot(melted_cormat, aes(Var2, Var1, fill = logfdr)) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colours = colors.martin) +
      geom_text(aes(Var2, Var1, label = round(fdr,3)), color = "black", size = 4) +
      xlab("") +
      ylab('') +
      ggtitle(title)
  }
}

#Cargar carpeta con otros scripts
#Hace un wilconson y un t test para la variable

R.utils::sourceDirectory("/mnt/beegfs/fheredia/Scripts_metilacion/Xec scripts/", modifiedOnly=F)

# Definir los grupos de compracion 
pdata <- targets
rownames(pdata) <- colnames(mSetSqFlt)

# We define the groups to compare
pdata_cov_DG <- pdata[,c(3,5,6,8,12,13,15,17,18)] #se eligen las columnas de las variables que se quieren testar
CCC <-  rownames(pdata_cov_DG)[pdata_cov_DG$Type == "Kuschnir"]
HD <-  rownames(pdata_cov_DG)[pdata_cov_DG$Type == "Healthy"]
K1<-  rownames(pdata_cov_DG)[pdata_cov_DG$Status == "Kuschnir1"]
K2<-  rownames(pdata_cov_DG)[pdata_cov_DG$Status == "Kuschnir2"]
K3<-  rownames(pdata_cov_DG)[pdata_cov_DG$Status == "Kuschnir3"]
K2K3<-  rownames(pdata_cov_DG)[pdata_cov_DG$Status3 == "K2K3"]

#colnames(mVals) #que coincidan con el nombre de las muestras
#memory.size(25000)

#Calcular covariables para cada comparacion 
#CCC vs HD
M_CCC_HD <- as.data.frame(mVals[,colnames(mVals) %in% c(CCC, HD)])
data_CCC_HD <- pdata_cov_DG[rownames(pdata_cov_DG) %in% c(CCC, HD),] 
matrix_CCC_HD <- pval_var(input = M_CCC_HD, pdata =data_CCC_HD, n_pca = 3, method = "fdr")
var_heatmap(input = matrix_CCC_HD, title = "CCC vs HD", pval = T)

#K2K3 vs HD
M_K2K3_HD <- as.data.frame(mVals[,colnames(mVals) %in% c(K2K3, HD)])
data_K2K3_HD <- pdata_cov_DG[rownames(pdata_cov_DG) %in% c(K2K3, HD),] 
matrix_K2K3_HD <- pval_var(input = M_K2K3_HD, pdata =data_K2K3_HD, n_pca = 3, method = "fdr")
var_heatmap(input = matrix_K2K3_HD, title = "K2K3 vs HD", pval = T)

#K1 vs HD
M_K1_HD <- as.data.frame(mVals[,colnames(mVals) %in% c(K1, HD)])
data_K1_HD <- pdata_cov_DG[rownames(pdata_cov_DG) %in% c(K1, HD),] 
matrix_K1_HD <- pval_var(input = M_K1_HD, pdata =data_K1_HD, n_pca = 3, method = "fdr")
var_heatmap(input = matrix_K1_HD, title = "K1 vs HD", pval = T)

#K2 vs HD
M_K2_HD <- as.data.frame(mVals[,colnames(mVals) %in% c(K2, HD)])
data_K2_HD <- pdata_cov_DG[rownames(pdata_cov_DG) %in% c(K2, HD),] 
matrix_K2_HD <- pval_var(input = M_K2_HD, pdata =data_K2_HD, n_pca = 3, method = "fdr")
var_heatmap(input = matrix_K2_HD, title = "K2 vs HD", pval = T)

#K3 vs HD
M_K3_HD <- as.data.frame(mVals[,colnames(mVals) %in% c(K3, HD)])
data_K3_HD <- pdata_cov_DG[rownames(pdata_cov_DG) %in% c(K3, HD),] 
matrix_K3_HD <- pval_var(input = M_K3_HD, pdata =data_K3_HD, n_pca = 3, method = "fdr")
var_heatmap(input = matrix_K3_HD, title = "K3 vs HD", pval = T)


#######################################################
###CALCULO DE SITIOS DIFERENCIALMENTE METILADOS (DMPs)  
#######################################################

#KUSCHNIR vs HEALTHY 
#design matrix
design<- model.matrix(~0+Status+predictedSex+Age, data=targets)
colnames(design) <- c("H", "K1", "K2", "K3", "M","Age")
design
#ajustar a modelo lineal
fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(K1_vs_H=K1-H,
                            K2_vs_H=K2-H,
                            K3_vs_H=K3-H,
                            levels=design)
contMatrix
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

 ### Kuschnir1 vs Healthy ###
DMPs_K1vsH <- topTable(fit2, num=Inf, coef="K1_vs_H")
head(DMPs_K1vsH)
#pasar datos a xlsx
#write.xlsx(DMPs_K1vsH, file = "DMPs_K1vsH.xlsx", rowNames=TRUE)

 ### Kuschnir2 vs Healthy ###
DMPs_K2vsH <- topTable(fit2, num=Inf, coef="K2_vs_H")
head(DMPs_K2vsH)
#pasar datos a xlsx
#write.xlsx(DMPs_K2vsH, file = "DMPs_K2vsH.xlsx", rowNames=TRUE)

 ### Kuschnir3 vs Healthy ###
DMPs_K3vsH <- topTable(fit2, num=Inf, coef="K3_vs_H")
head(DMPs_K3vsH)
#pasar datos a xlsx
#write.xlsx(DMPs_K3vsH, file = "DMPs_K3vsH.xlsx", rowNames=TRUE)


  #K1 K2 Y K3 JUNTOS (CCC) VS HEALTHY 
design_2<- model.matrix(~0+Type+predictedSex+Age, data=targets)
colnames(design_2) <- c("H", "CCC", "M","Age")
design_2
fit_2 <- lmFit(mVals, design_2)

contMatrix_2 <- makeContrasts(CCC_vs_H=CCC-H,
                              levels=design_2)
contMatrix_2
fit2_2 <- contrasts.fit(fit_2, contMatrix_2)
fit2_2 <- eBayes(fit2_2)
summary(decideTests(fit2_2))

DMPs_CCCvsH <- topTable(fit2_2, num=Inf, coef="CCC_vs_H")
head(DMPs_CCCvsH)
#pasar datos a xlsx
#write.xlsx(DMPs_CCCvsH, file = "DMPs_CCCvsH.xlsx", rowNames=TRUE)


  #K2K3 VS HEALTHY
#Para hacerlo creo nueva columna en targets
targets$Status3 <- ifelse(targets$Status %in% c("Kuschnir1", "Healthy"), targets$Status, "K2K3")

design_3<- model.matrix(~0+Status3+predictedSex+Age, data=targets)
colnames(design_3) <- c("H", "K2K3","K1", "M","Age")
design_3
fit_3 <- lmFit(mVals, design_3)

contMatrix_3 <- makeContrasts(K2K3_vs_H=K2K3-H,
                              levels=design_3)
contMatrix_3
fit2_3 <- contrasts.fit(fit_3, contMatrix_3)
fit2_3 <- eBayes(fit2_3)
summary(decideTests(fit2_3))
DMPs_K2K3vsH <- topTable(fit2_3, num=Inf, coef="K2K3_vs_H")
head(DMPs_K2K3vsH)
#pasar datos a xlsx
#write.xlsx(DMPs_K2K3vsH, file = "DMPs_K2K3vsH.xlsx", rowNames=TRUE)


#############################################
### DIFFERENTIAL VARIABILITY POSITIONS (DVPs)
#############################################
#https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html -- punto 3.2

#en este caso se hace un fit a la varianza de los datos (varFit)
## K1, K2 y K3 vs HD
fitvar <- varFit(mVals, design = design)
fitvar2 <- contrasts.varFit(fitvar, contMatrix)
summary(decideTests(fitvar2))
DVPs_K1vsH <- topVar(fitvar2, coef="K1_vs_H", num=300)
DVPs_K2vsH <- topVar(fitvar2, coef="K2_vs_H", num=300)
DVPs_K3vsH <- topVar(fitvar2, coef="K3_vs_H", num=300)

## CCC vs HD
fitvar_2 <- varFit(mVals, design = design_2)
fitvar2_2 <- contrasts.varFit(fitvar_2, contMatrix_2)
summary(decideTests(fitvar2_2))
DVPs_CCCvsH <- topVar(fitvar2_2, coef="CCC_vs_H", num=300)

## K2K3 vs HD
fitvar_3 <- varFit(mVals, design = design_3)
fitvar2_3 <- contrasts.varFit(fitvar_3, contMatrix_3)
summary(decideTests(fitvar2_3))
DVPs_K2K3vsH <- topVar(fitvar2_3, coef="K2K3_vs_H", num=300)


### REPRESENTACION DE LAS TOP 4 DVPs
setwd("/mnt/beegfs/fheredia/Scripts_metilacion/")
pdf(file=paste("DVPs", "_", format(Sys.Date(), "%d.%m.%Y"), ".pdf", sep = ""), width=9, height=8, paper="a4r")

# K1 vs HD
par(mfrow=c(2,2), oma=c(1,1,2,1))
sapply(rownames(DVPs_K1vsH)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Status, 
          ylab = "Beta values")
})
mtext("Top 4 DVPs K1 vs HD", outer = TRUE, cex = 1.2)
# K2 vs HD
par(mfrow=c(2,2), oma=c(1,1,2,1))
sapply(rownames(DVPs_K2vsH)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Status, 
          ylab = "Beta values")
})
mtext("Top 4 DVPs K2 vs HD", outer = TRUE, cex = 1.2)
# K3 vs HD
par(mfrow=c(2,2), oma=c(1,1,2,1))
sapply(rownames(DVPs_K3vsH)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Status, 
          ylab = "Beta values")
})
mtext("Top 4 DVPs K3 vs HD", outer = TRUE, cex = 1.2)
# CCC vs HD
par(mfrow=c(2,2), oma=c(1,1,2,1))
sapply(rownames(DVPs_CCCvsH)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Type, 
          ylab = "Beta values")
})
mtext("Top 4 DVPs CCC vs HD", outer = TRUE, cex = 1.2)
# K2K3 vs HD
par(mfrow=c(2,2), oma=c(1,1,2,1))
sapply(rownames(DVPs_K2K3vsH)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Status3, 
          ylab = "Beta values")
})
mtext("Top 4 DVPs K2K3 vs HD", outer = TRUE, cex = 1.2)

dev.off()


#to be continued...
