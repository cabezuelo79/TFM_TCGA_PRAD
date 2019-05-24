##############################################################
####### Analisis de expresion genica mediante RNA-Seq#########        
############## ESTUDIO PAREADO ###############################
###### Analisis de expresion diferencial con el paquete edgeR#
##############################################################

# Establecimiento de directorio de trabajo, donde se encuentra el archivo "PAREADO.txt"
setwd('C:/Users/JoseMaria/OneDrive - MERIDIEM SEEDS S.L/Documentos/Master/TFM/RNASeq_PRAD')

# Cargo el paquete necesario para la ejecución del analisis
library(edgeR)

# Obtengo la información del archivo "PAREADO.txt" 
rawdataPAR <- read.delim("PAREADO.txt", check.names=FALSE, stringsAsFactors=FALSE)

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes, el resto las muestras
yPAR <- DGEList(counts=rawdataPAR[,3:106], genes=rawdataPAR[,1:2])

# Llevo a cabo la normalizacion de las muestras
yPAR <- calcNormFactors(yPAR)

# Cargo los vectores 'Tissue' y 'Patient' de forma manual
TissuePAR <- factor(c("N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T"))
PatientPAR <- factor(c("CH-5761","CH-5761","CH-5767","CH-5767","CH-5768","CH-5768","CH-5769","CH-5769","EJ-7115","EJ-7115","EJ-7123","EJ-7123","EJ-7125","EJ-7125","EJ-7314","EJ-7314","EJ-7315","EJ-7315","EJ-7317","EJ-7317","EJ-7321","EJ-7321","EJ-7327","EJ-7327","EJ-7328","EJ-7328","EJ-7330","EJ-7330","EJ-7331","EJ-7331","EJ-7781","EJ-7781","EJ-7782","EJ-7782","EJ-7783","EJ-7783","EJ-7784","EJ-7784","EJ-7785","EJ-7785","EJ-7786","EJ-7786","EJ-7789","EJ-7789","EJ-7792","EJ-7792","EJ-7793","EJ-7793","EJ-7794","EJ-7794","EJ-7797","EJ-7797","EJ-A8FO","EJ-A8FO","G9-6333","G9-6333","G9-6342","G9-6342","G9-6348","G9-6348","G9-6351","G9-6351","G9-6356","G9-6356","G9-6362","G9-6362","G9-6363","G9-6363","G9-6365","G9-6365","G9-6384","G9-6384","G9-6496","G9-6496","G9-6499","G9-6499","HC-7211","HC-7211","HC-7737","HC-7737","HC-7738","HC-7738","HC-7740","HC-7740","HC-7742","HC-7742","HC-7745","HC-7745","HC-7747","HC-7747","HC-7752","HC-7752","HC-7819","HC-7819","HC-8258","HC-8258","HC-8259","HC-8259","HC-8260","HC-8260","HC-8262","HC-8262","J4-A83J","J4-A83J"))

# Construyo un data frame con la informacion
data.frame(Sample=colnames(yPAR),PatientPAR,TissuePAR)

# Construyo la matriz de diseño en base a los dos factores ('Tissue' y 'Patient')
designPAR <- model.matrix(~PatientPAR+TissuePAR)
rownames(designPAR) <- colnames(yPAR)

# Calculo la dispersion estimada y visualizo la dispersion comun
yPAR <- estimateDisp(yPAR, designPAR, robust=TRUE) 
yPAR$common.dispersion

# Visualizo la grafica del coeficiente de variacion biologica (BCV)
plotBCV(yPAR)

# Likelihood-ratio test
fitPAR <- glmFit(yPAR, designPAR)
lrtPAR <- glmLRT(fitPAR)

# Visualizo un plot multidimensional con los genes sub-expresados y sobre-expresados
plotMD(lrtPAR)

# Genero un elemento 'topTags' con el resultado obtenido (n=Inf para que incluya todos los resultados)
top_PAR <- topTags(lrtPAR,n=Inf)

# visualizo el resultado
top_PAR

# visualizo el resumen de genes sobre-expresados, genes sub-expresados y genes no significativos
summary(decideTests(lrtPAR))

##############################################################################################
############## COMPROBACION DE RESULTADOS OBTENIDOS MANUALMENTE ##############################
##############################################################################################

# Construyo el data frame 'TagetsPAR' con los counts y el tejido correspondiente
targetsPAR<-data.frame(filename=colnames(rawdataPAR[,3:106]),type=TissuePAR)

# Visualizo 'targetsPAR'
head(targetsPAR)

# Cambio 'filename' de factor a caracter
targetsPAR$filename<-as.character(targetsPAR$filename)

# Genero la media de counts para el primero de los genes diferencialmente expresados y tejido normal
controlPAR <- mean(as.numeric(yPAR$counts[yPAR$genes$`SYMBOL`=="TRIP6",targetsPAR[targetsPAR$type=="N","filename"]])) 

# Visualizo el resultado
controlPAR

# Genero la media de counts para el primero de los genes diferencialmente expresados y tejido tumoral
samplePAR <- mean(as.numeric(yPAR$counts[yPAR$genes$`SYMBOL`=="TRIP6",targetsPAR[targetsPAR$type=="T","filename"]])) 

# Visualizo el resultado
samplePAR

# Realizo el calculo de logFC manualmente
log2(samplePAR/controlPAR)



##############################################################################################
############## ANALISIS DE ENRIQUECIMIENTO FUNCIONAL GOSEQ      ##############################
##############################################################################################

# Cargo los paquetes necesarios para el analisis
library(goseq)
library(GO.db)

# Genero un vector con los genes diferencialmente expresados (FDR < 0.05)
topPAR_DE <- top_PAR$table[top_PAR$table$FDR < 0.05,"SYMBOL"]

# Genero un vector 'universal' que contiene tanto los genes DE como los no DE
universePAR <- top_PAR$table$`SYMBOL`

# Construyo un vector apropiado para el analisis con 'goseq'
mytablePAR <- as.integer(unique(universePAR)%in%topPAR_DE)
names(mytablePAR) <- unique(universePAR)

# Veo la cantidad de genes DE (1) y no DE (0) y compruebo que este correcto
table(mytablePAR)
head(mytablePAR)

# Probability Weighting Function (pwf) utilizando 'geneSymbol' y genoma 'hg19'
pwf_PAR=nullp(mytablePAR,"hg19","geneSymbol")

# Compruebo el resultado
head(pwf_PAR)

# Realizo el analisis GO mediante la aproximación de Wallenius
GO.wall_PAR=goseq(pwf_PAR,"hg19","geneSymbol")

# Visualizo el resultado
GO.wall_PAR

# Llevo a cabo el enriquecimiento
enriched.GO_PAR <- GO.wall_PAR$category[p.adjust(GO.wall_PAR$over_represented_pvalue, method="BH")<.05]

# Visualizo el resultado
for(go in enriched.GO_PAR[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

##### ANALISIS DE RUTAS KEGG #########

kegg_DE_PAR <- goseq(pwf_PAR,'hg19','geneSymbol',test.cats="KEGG")

# Visualizo el resultado
kegg_DE_PAR
