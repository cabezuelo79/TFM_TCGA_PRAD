##############################################################
####### Analisis de expresion genica mediante RNA-Seq#########        
############## ESTUDIO GENERAL ###############################
###### Analisis de expresion diferencial con el paquete edgeR#
##############################################################

# Establecimiento de directorio de trabajo, donde se encuentra el archivo "GENERAL.txt"
dirJM<-"C:/Users/JoseMaria/OneDrive - MERIDIEM SEEDS S.L/Documentos/Master/TFM/RNASeq_PRAD/GitHub/TFM_TCGA_PRAD"
dirEdu<-"/Users/eandres/Proyectos/Eduardo_Andres/TFM_Cabezuelo/Cabezuelo/TFM_TCGA_PRAD"
setwd(dirEdu)
setwd(dirJM)

mapPathwayToName <- function(organism) {
  KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
  pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, organism, sep="")
  
  pathway_id_name <- data.frame()
  cont<-0
  for (line in readLines(pathway_list_REST_url)) {
    cont<-cont+1
    tmp <- strsplit(line, "\t")[[1]]
    pathway_id <- strsplit(tmp[1], organism)[[1]][2]
    pathway_name <- tmp[2]
    pathway_name <- strsplit(pathway_name, "\\s+-\\s+")[[1]][1]
    pathway_id_name[cont, 1] = pathway_id
    pathway_id_name[cont, 2] = pathway_name
    
  }
  names(pathway_id_name) <- c("path","pathway_name")
  pathway_id_name
}

# Cargo el paquete necesario para la ejecución del analisis
library(edgeR)

# Obtengo la información del archivo "GENERAL.txt" 
rawdata <- read.delim("GENERAL.txt", check.names=FALSE, stringsAsFactors=FALSE)

#incluyo el nombre de los genes como rowname
rownames(rawdata)<-rawdata$SYMBOL

# Cargo el vector 'Tissue' de forma manual
Tissue <- c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","N","T","N","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","N","T","N","T","T","T","N","T","N","T","N","T","T","N","T","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","T","N","T","T","N","T","N","T","N","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","N","T","T","T","N","T","N","T","T","T","N","T","T","N","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","N","T","T","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","N","T","N","T","N","T","N","T","T","N","T","N","T","T","T","T","N","T","T","T","N","T","T","T","T","T","T","T","N","T","N","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T")

#Obtengo el tamaño de cada gen para normalizar teniendo en ceunta del tamaño de los genes
gene.length<-read.table("his-Size.tab",header=T)
idx<-match(rawdata$SYMBOL,gene.length$Gene)
results_counts<-gene.length[idx,]
results_counts[is.na(results_counts$Length),"Length"]<-0
nrow(results_counts)

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes (NO!! eso era el ID segun Refseq), el resto las muestras
# como grupo selecciono tipo de tejido 'Tissue'
yET <- DGEList(counts=rawdata[,3:551], genes=results_counts, group=Tissue)

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
head(top$table)

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

head(yET$counts)

# Genero la media de counts para el primero de los genes diferencialmente expresados y tejido normal
control <- mean(as.numeric(yET$counts["SERPINA5",targets[targets$type=="N","filename"]]))

# Visualizo el resultado
control

# Genero la media de counts para el primero de los genes diferencialmente expresados y tejido tumoral
sample <- mean(as.numeric(yET$counts["SERPINA5",targets[targets$type=="T","filename"]]))

# Visualizo el resultado
sample

# Realizo el calculo de logFC manualmente
log2(sample/control)
top$table["SERPINA5","logFC"]

##############################################################################################
############## ANALISIS DE ENRIQUECIMIENTO FUNCIONAL GOSEQ      ##############################
##############################################################################################

# Cargo los paquetes necesarios para el analisis
library(goseq)
library(GO.db)

# Genero un vector con los genes diferencialmente expresados (FDR < 0.05)
topET_DE <- as.vector(na.omit(top$table[top$table$FDR < 0.05,"Gene"]))

# Genero un vector 'universal' que contiene tanto los genes DE como los no DE
universe <- as.vector(unique(na.omit(top$table$Gene)))

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
GO.wall=goseq(pwf,"hg19","geneSymbol", test.cats = c("GO:BP","GO:CC","GO:MF"))

# Separo los terminos GO y creo FDRdown y FDRover
GO.BP <- GO.wall[GO.wall$ontology=="BP",]
GO.BP$FDRunder <- p.adjust(GO.BP$under_represented_pvalue, n=nrow(GO.BP))
GO.BP$FDRover <- p.adjust(GO.BP$over_represented_pvalue, n=nrow(GO.BP))
write.table(GO.BP, file = 'GO.biological.process.tsv', sep = "\t", row.names = FALSE)

GO.CC <- GO.wall[GO.wall$ontology=="CC",]
GO.CC$FDRunder <- p.adjust(GO.CC$under_represented_pvalue, n=nrow(GO.CC))
GO.CC$FDRover <- p.adjust(GO.CC$over_represented_pvalue, n=nrow(GO.CC))
write.table(GO.CC, file = 'GO.cellular.component.tsv', sep = "\t", row.names = FALSE)

GO.MF <- GO.wall[GO.wall$ontology=="MF",]
GO.MF$FDRunder <- p.adjust(GO.MF$under_represented_pvalue, n=nrow(GO.MF))
GO.MF$FDRover <- p.adjust(GO.MF$over_represented_pvalue, n=nrow(GO.MF))
write.table(GO.MF, file = 'GO.molecular.function.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP <- GO.BP[GO.BP$FDRover<=0.05,]
enriched.over.GO.CC <- GO.CC[GO.CC$FDRover<=0.05,]
enriched.over.GO.MF <- GO.MF[GO.MF$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP <- GO.BP[GO.BP$FDRunder<=0.05,]
write.table(enriched.under.GO.BP, file = 'GO.biological.process_under.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC <- GO.CC[GO.CC$FDRunder<=0.05,]
write.table(enriched.under.GO.CC, file = 'GO.cellular.component_under.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF <- GO.MF[GO.MF$FDRunder<=0.05,]
write.table(enriched.under.GO.MF, file = 'GO.molecular.function_under.tsv', sep = "\t", row.names = FALSE)




##### ANALISIS DE RUTAS KEGG #########

kegg_DE <- goseq(pwf,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE$FDRunder <- p.adjust(kegg_DE$under_represented_pvalue, n=nrow(kegg_DE))
kegg_DE$FDRover <- p.adjust(kegg_DE$over_represented_pvalue, n=nrow(kegg_DE))

kegg_path<-mapPathwayToName("hsa")
idx<-match(kegg_DE$category,kegg_path$path)
pathway_name<-kegg_path[idx,2]
data<-cbind(kegg_DE,pathway_name)
# Visualizo el resultado
head(data)

write.table(data, file = 'KEGG.PATHWAYS.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG <- data[data$FDRover<=0.05,]

# Genes under-expressed
under.KEGG <- data[data$FDRunder<=0.05,]
write.table(under.KEGG, file = 'KEGG.PATHWAYS_under.tsv', sep = "\t", row.names = FALSE)


