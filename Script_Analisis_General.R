##############################################################
####### Analisis de expresion genica mediante RNA-Seq#########        
############## ESTUDIO GENERAL ###############################
###### Analisis de expresion diferencial con el paquete edgeR#
##############################################################

# Establecimiento de directorio de trabajo, donde se encuentra el archivo "GENERAL.txt"
setwd('C:/Users/JoseMaria/OneDrive - MERIDIEM SEEDS S.L/Documentos/Master/TFM/RNASeq_PRAD')

# Cargo el paquete necesario para la ejecución del analisis
library(edgeR)

# Obtengo la información del archivo "GENERAL.txt" 
rawdata <- read.delim("GENERAL.txt", check.names=FALSE, stringsAsFactors=FALSE)

# Cargo el vector 'Tissue' de forma manual
Tissue <- c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","N","T","N","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","N","T","N","T","T","T","N","T","N","T","N","T","T","N","T","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","T","N","T","T","N","T","N","T","N","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","N","T","T","T","N","T","N","T","T","T","N","T","T","N","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","N","T","T","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","N","T","N","T","N","T","N","T","T","N","T","N","T","T","T","T","N","T","T","T","N","T","T","T","T","T","T","T","N","T","N","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T"
)

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes, el resto las muestras
# como grupo selecciono tipo de tejido 'Tissue'
yET <- DGEList(counts=rawdata[,3:551], genes=rawdata[,1:2], group=Tissue)

# Llevo a cabo la normalizacion de las muestras
dgenorm<-calcNormFactors(yET)

# Realizo la estimacion de la dispersión y BCV
dgenorm <- estimateCommonDisp(dgenorm,verbose = T)

# Realizo la dispersion 'tagwise' despues del calculo de la dispersion comun
dgenorm <- estimateTagwiseDisp(dgenorm)

# Llevo a cabo el analisis de expresion diferencial mediante 'exactTest' con el par 'N' tejido normal y 'T' tejido tumoral
et <- exactTest(dgenorm,pair=c("N","T"))

# Genero una tabla 'topTags' con la información obtenida
top<-topTags(et, n=nrow(et), adjust.method="BH",sort.by="PValue")

# visualizo el resultado
top

# visualizo el resumen de genes sobre-expresados, genes sub-expresados y genes no significativos
summary(decideTests(et))


##############################################################################################
############## COMPROBACION DE RESULTADOS OBTENIDOS MANUALMENTE ##############################
##############################################################################################

# Construyo el data frame 'Tagets' con los counts y el tejido correspondiente
targets<-data.frame(filename=colnames(rawdata[,3:551]),type=Tissue)

# Visualizo 'targets'
head(targets)

# Cambio 'filename' de factor a caracter
targets$filename<-as.character(targets$filename)

# Genero la media de counts para el primero de los genes diferencialmente expresados y tejido normal
control <- mean(as.numeric(yET$counts[yET$genes$`SYMBOL`=="SERPINA5",targets[targets$type=="N","filename"]]))

# Visualizo el resultado
control

# Genero la media de counts para el primero de los genes diferencialmente expresados y tejido tumoral
sample <- mean(as.numeric(yET$counts[yET$genes$`SYMBOL`=="SERPINA5",targets[targets$type=="T","filename"]]))

# Visualizo el resultado
sample

# Realizo el calculo de logFC manualmente
log2(sample/control)


##############################################################################################
############## ANALISIS DE ENRIQUECIMIENTO FUNCIONAL GOSEQ      ##############################
##############################################################################################

# Cargo los paquetes necesarios para el analisis
library(goseq)
library(GO.db)

# Genero un vector con los genes diferencialmente expresados (FDR < 0.05)
topET_DE <- top$table[top$table$FDR < 0.05,"SYMBOL"]

# Genero un vector 'universal' que contiene tanto los genes DE como los no DE
universe <- top$table$`SYMBOL`

# Construyo un vector apropiado para el analisis con 'goseq'
mytable <- as.integer(unique(universe)%in%topET_DE)
names(mytable) <- unique(universe)

# Veo la cantidad de genes DE (1) y no DE (0) y compruebo que este correcto
table(mytable)
head(mytable)

# Probability Weighting Function (pwf) utilizando 'geneSymbol' y genoma 'hg19'
pwf=nullp(mytable,"hg19","geneSymbol")

# Compruebo el resultado
head(pwf)

# Realizo el analisis GO mediante la aproximación de Wallenius
GO.wall=goseq(pwf,"hg19","geneSymbol")

# Visualizo el resultado
GO.wall

# Llevo a cabo el enriquecimiento
enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]

# Visualizo el resultado
for(go in enriched.GO[1:10]){
      print(GOTERM[[go]])
      cat("--------------------------------------\n")
      }

##### ANALISIS DE RUTAS KEGG #########

kegg_DE <- goseq(pwf,'hg19','geneSymbol',test.cats="KEGG")

# Visualizo el resultado
kegg_DE



