##############################################################
####### Analisis de expresion genica mediante RNA-Seq#########        
##############       GLEASON   ###############################
###### Analisis de expresion diferencial con el paquete edgeR#
##############################################################

# Establecimiento de directorio de trabajo, donde se encuentra el archivo "GENERAL.txt"
setwd('C:/Users/JoseMaria/OneDrive - MERIDIEM SEEDS S.L/Documentos/Master/TFM/RNASeq_PRAD/GitHub/TFM_TCGA_PRAD')

# Cargo el paquete necesario para la ejecución del analisis
library(edgeR)

# Obtengo la información del archivo "GENERAL.txt" 
rawdata <- read.delim("GENERAL.txt", check.names=FALSE, stringsAsFactors=FALSE)

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes, el resto las muestras
yG <- DGEList(counts=rawdata[,3:106], genes=rawdata[,1:2])

# Cargo los vectores 'Sample' y 'Gleason' de forma manual
Gleason <- factor(c("G1","G1","G2","G1","G2","G2","G2","G2","G1","G1","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G2","G2","G2","G2","G2","G0","G2","G0","G1","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G0","G2","G0","G2","G2","G2","G0","G2","G0","G2","G0","G2","G2","G0","G1","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G0","G2","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G1","G2","G2","G0","G2","G2","G2","G2","G0","G1","G2","G1","G0","G2","G0","G2","G2","G2","G0","G2","G2","G0","G2","G0","G2","G2","G0","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G0","G2","G2","G2","G0","G2","G2","G0","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G1","G1","G2","G2","G2","G2","G1","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G0","G2","G0","G2","G1","G2","G2","G0","G2","G2","G2","G0","G2","G2","G2","G1","G2","G2","G2","G0","G1","G0","G1","G0","G2","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G0","G2","G1","G2","G2","G2","G1","G2","G1","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G1","G2","G2","G1","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2"
))

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes, el resto las muestras
# como grupo selecciono tipo de tejido 'Gleason'
yG <- DGEList(counts=rawdata[,3:551], genes=rawdata[,1:2], group=Gleason)

# Llevo a cabo la normalizacion de las muestras
yG <- calcNormFactors(yG)

# Examino las muestras mediante grafico MDS
plotMDS(yG)

#Diseño la matriz para glm
designG <- model.matrix(~0+Gleason, data=yG$samples)
colnames(designG) <- levels(yG$samples$group)

# visualizo la matriz
designG

# Calculo las dispersiones
yG <- estimateGLMCommonDisp(yG, designG)
yG <- estimateGLMTrendedDisp(yG, designG)
yG <- estimateGLMTagwiseDisp(yG, designG)

# visualizo la dispersion en un grafico BCV
plotBCV(yG)

# Quasi-likelihood test
fitG <- glmQLFit(yG, designG)

#CONTRASTES
my.contrasts <- makeContrasts(G1vsG0=G1-G0, G2vsG1=G2-G1, G2vsG0=G2-G0, levels=designG)
qlf.G1vsG0 <- glmQLFTest(fitG, contrast=my.contrasts[,"G1vsG0"])
top_G1vsG0 <- topTags(qlf.G1vsG0, n=Inf)
qlf.G2vsG1 <- glmQLFTest(fitG, contrast=my.contrasts[,"G2vsG1"])
top_G2vsG1 <- topTags(qlf.G2vsG1, n=Inf)
qlf.G2vsG0 <- glmQLFTest(fitG, contrast=my.contrasts[,"G2vsG0"])
top_G2vsG0 <- topTags(qlf.G2vsG0, n=Inf)

# visualizo el resumen de genes sobre-expresados, genes sub-expresados y genes no significativos
summary(decideTests(qlf.G1vsG0))
summary(decideTests(qlf.G2vsG1))
summary(decideTests(qlf.G2vsG0))

##############################################################################################
############## ANALISIS DE ENRIQUECIMIENTO FUNCIONAL GOSEQ      ##############################
##############################################################################################

# Cargo los paquetes necesarios para el analisis
library(goseq)
library(GO.db)

# Genero vectores con los genes diferencialmente expresados (FDR < 0.05)
top_G1vsG0_filter <- top_G1vsG0$table[top_G1vsG0$table$FDR < 0.05,"SYMBOL"]
top_G2vsG1_filter <- top_G2vsG1$table[top_G2vsG1$table$FDR < 0.05,"SYMBOL"]
top_G2vsG0_filter <- top_G2vsG0$table[top_G2vsG0$table$FDR < 0.05,"SYMBOL"]

# Genero un vector 'universal' que contiene tanto los genes DE como los no DE
universeG1vsG0 <- top_G1vsG0$table$`SYMBOL`
universeG2vsG1 <- top_G2vsG1$table$`SYMBOL`
universeG2vsG0 <- top_G2vsG0$table$`SYMBOL`

# Construyo vectores apropiados para el analisis con 'goseq'
mytableG1vsG0 <- as.integer(unique(universeG1vsG0)%in%top_G1vsG0_filter)
names(mytableG1vsG0) <- unique(universeG1vsG0)
mytableG2vsG1 <- as.integer(unique(universeG2vsG1)%in%top_G2vsG1_filter)
names(mytableG2vsG1) <- unique(universeG2vsG1)
mytableG2vsG0 <- as.integer(unique(universeG2vsG0)%in%top_G2vsG0_filter)
names(mytableG2vsG0) <- unique(universeG2vsG0)

# Veo la cantidad de genes DE (1) y no DE (0) y compruebo que este correcto
table(mytableG1vsG0)
head(mytableG1vsG0)
table(mytableG2vsG1)
head(mytableG2vsG1)
table(mytableG2vsG0)
head(mytableG2vsG0)


# Probability Weighting Function (pwf) utilizando 'geneSymbol' y genoma 'hg19'
pwf_G1vsG0=nullp(mytableG1vsG0,"hg19","geneSymbol")
pwf_G2vsG1=nullp(mytableG2vsG1,"hg19","geneSymbol")
pwf_G2vsG0=nullp(mytableG2vsG0,"hg19","geneSymbol")

# Compruebo el resultado
head(pwf_G1vsG0)
head(pwf_G2vsG1)
head(pwf_G2vsG0)


# Realizo el analisis GO mediante la aproximación de Wallenius
GO.wall_G1vsG0=goseq(pwf_G1vsG0,"hg19","geneSymbol")
GO.wall_G2vsG1=goseq(pwf_G2vsG1,"hg19","geneSymbol")
GO.wall_G2vsG0=goseq(pwf_G2vsG0,"hg19","geneSymbol")

# Visualizo el resultado
GO.wall_G1vsG0
GO.wall_G2vsG1
GO.wall_G2vsG0

# Llevo a cabo el enriquecimiento
enriched.GO_G1vsG0 <- GO.wall_G1vsG0$category[p.adjust(GO.wall_G1vsG0$over_represented_pvalue, method="BH")<.05]
enriched.GO_G2vsG1 <- GO.wall_G2vsG1$category[p.adjust(GO.wall_G2vsG1$over_represented_pvalue, method="BH")<.05]
enriched.GO_G2vsG0 <- GO.wall_G2vsG0$category[p.adjust(GO.wall_G2vsG0$over_represented_pvalue, method="BH")<.05]

# Visualizo el resultado
for(go in enriched.GO_G1vsG0[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

for(go in enriched.GO_G2vsG1[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

for(go in enriched.GO_G2vsG0[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

##### ANALISIS DE RUTAS KEGG #########

kegg_DE_G1vsG0 <- goseq(pwf_G1vsG0,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_G2vsG1 <- goseq(pwf_G2vsG1,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_G2vsG0 <- goseq(pwf_G2vsG0,'hg19','geneSymbol',test.cats="KEGG")

# Visualizo el resultado
kegg_DE_G1vsG0
kegg_DE_G2vsG1
kegg_DE_G2vsG0