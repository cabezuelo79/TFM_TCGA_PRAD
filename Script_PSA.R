##############################################################
####### Analisis de expresion genica mediante RNA-Seq#########        
##############       PSA       ###############################
###### Analisis de expresion diferencial con el paquete edgeR#
##############################################################

# Establecimiento de directorio de trabajo, donde se encuentra el archivo "psa.txt"
setwd('C:/Users/Usuario/OneDrive - MERIDIEM SEEDS S.L/Documentos/Master/TFM/RNASeq_PRAD/GitHub/TFM_TCGA_PRAD')

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

# Obtengo la información del archivo "psa.txt" 
rawdata_psa <- read.delim("psa.txt", check.names=FALSE, stringsAsFactors=FALSE)

#incluyo el nombre de los genes como rowname
rownames(rawdata_psa)<-rawdata_psa$SYMBOL

# Cargo el vector 'Gleason' de forma manual
psa <- factor(c("PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA1","PSA0","PSA1","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA2","PSA0","PSA1","PSA1","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA1","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA1","PSA1","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA2","PSA0","PSA1","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0"))

#Obtengo el tamaño de cada gen para normalizar teniendo en ceunta del tamaño de los genes
gene.length_psa <-read.table("his-Size.tab",header=T)
idx_psa <-match(rawdata_psa$SYMBOL,gene.length_psa$Gene)
results_counts_psa <-gene.length_psa[idx_psa,]
results_counts_psa[is.na(results_counts_psa$Length),"Length"]<-0
nrow(results_counts_psa)

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes, el resto las muestras
# como grupo selecciono tipo de tejido 'Gleason'
yPSA <- DGEList(counts=rawdata_psa[,3:442], genes=results_counts_psa, group=psa)

# Llevo a cabo la normalizacion de las muestras
yPSA <- calcNormFactors(yPSA)

# Examino las muestras mediante grafico MDS
plotMDS(yPSA)

#Diseño la matriz para glm
designPSA <- model.matrix(~0+psa, data=yPSA$samples)
colnames(designPSA) <- levels(yPSA$samples$group)

# visualizo la matriz
designPSA

# Calculo las dispersiones
yPSA <- estimateGLMCommonDisp(yPSA, designPSA)
yPSA <- estimateGLMTrendedDisp(yPSA, designPSA)
yPSA <- estimateGLMTagwiseDisp(yPSA, designPSA)

# visualizo la dispersion en un grafico BCV
plotBCV(yPSA)

# Quasi-likelihood test
fitPSA <- glmQLFit(yPSA, designPSA)

#CONTRASTES
my.contrasts_PSA <- makeContrasts(PSA1vsPSA0=PSA1-PSA0, PSA2vsPSA1=PSA2-PSA1, PSA2vsPSA0=PSA2-PSA0, levels=designPSA)
qlf.PSA1vsPSA0 <- glmQLFTest(fitPSA, contrast=my.contrasts_PSA[,"PSA1vsPSA0"])
top_PSA1vsPSA0 <- topTags(qlf.PSA1vsPSA0, n=Inf)
qlf.PSA2vsPSA1 <- glmQLFTest(fitPSA, contrast=my.contrasts_PSA[,"PSA2vsPSA1"])
top_PSA2vsPSA1 <- topTags(qlf.PSA2vsPSA1, n=Inf)
qlf.PSA2vsPSA0 <- glmQLFTest(fitPSA, contrast=my.contrasts_PSA[,"PSA2vsPSA0"])
top_PSA2vsPSA0 <- topTags(qlf.PSA2vsPSA0, n=Inf)

# visualizo el resumen de genes sobre-expresados, genes sub-expresados y genes no significativos
summary(decideTests(qlf.PSA1vsPSA0))
summary(decideTests(qlf.PSA2vsPSA1))
summary(decideTests(qlf.PSA2vsPSA0))

##############################################################################################
############## ANALISIS DE ENRIQUECIMIENTO FUNCIONAL GOSEQ      ##############################
##############################################################################################

# Cargo los paquetes necesarios para el analisis
library(goseq)
library(GO.db)

# Genero vectores con los genes diferencialmente expresados (FDR < 0.05)
top_PSA1vsPSA0_filter <- as.vector(na.omit(top_PSA1vsPSA0$table[top_PSA1vsPSA0$table$FDR < 0.05,"Gene"]))
top_PSA2vsPSA1_filter <- as.vector(na.omit(top_PSA2vsPSA1$table[top_PSA2vsPSA1$table$FDR < 0.05,"Gene"]))
top_PSA2vsPSA0_filter <- as.vector(na.omit(top_PSA2vsPSA0$table[top_PSA2vsPSA0$table$FDR < 0.05,"Gene"]))

# Genero vectores 'universal' que contienen tanto los genes DE como los no DE
universePSA1vsPSA0 <- as.vector(unique(na.omit(top_PSA1vsPSA0$table$Gene)))
universePSA2vsPSA1 <- as.vector(unique(na.omit(top_PSA2vsPSA1$table$Gene)))
universePSA2vsPSA0 <- as.vector(unique(na.omit(top_PSA2vsPSA0$table$Gene)))

# Construyo vectores apropiados para el analisis con 'goseq'
mytablePSA1vsPSA0 <- as.integer(unique(universePSA1vsPSA0)%in%top_PSA1vsPSA0_filter)
names(mytablePSA1vsPSA0) <- unique(universePSA1vsPSA0)
mytablePSA2vsPSA1 <- as.integer(unique(universePSA2vsPSA1)%in%top_PSA2vsPSA1_filter)
names(mytablePSA2vsPSA1) <- unique(universePSA2vsPSA1)
mytablePSA2vsPSA0 <- as.integer(unique(universePSA2vsPSA0)%in%top_PSA2vsPSA0_filter)
names(mytablePSA2vsPSA0) <- unique(universePSA2vsPSA0)

# Veo la cantidad de genes DE (1) y no DE (0) y compruebo que este correcto
table(mytablePSA1vsPSA0)
head(mytablePSA1vsPSA0)
table(mytablePSA2vsPSA1)
head(mytablePSA2vsPSA1)
table(mytablePSA2vsPSA0)
head(mytablePSA2vsPSA0)


# Probability Weighting Function (pwf) utilizando 'geneSymbol' y genoma 'hg19'
pwf_PSA1vsPSA0=nullp(mytablePSA1vsPSA0,"hg19","geneSymbol")
pwf_PSA2vsPSA1=nullp(mytablePSA2vsPSA1,"hg19","geneSymbol")
pwf_PSA2vsPSA0=nullp(mytablePSA2vsPSA0,"hg19","geneSymbol")

# Compruebo el resultado
head(pwf_PSA1vsPSA0)
head(pwf_PSA2vsPSA1)
head(pwf_PSA2vsPSA0)


# Realizo los analisis GO mediante la aproximación de Wallenius
GO.wall_PSA1vsPSA0=goseq(pwf_PSA1vsPSA0,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
GO.wall_PSA2vsPSA1=goseq(pwf_PSA2vsPSA1,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
GO.wall_PSA2vsPSA0=goseq(pwf_PSA2vsPSA0,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
kegg_DE_PSA1vsPSA0 <- goseq(pwf_PSA1vsPSA0,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_PSA2vsPSA1 <- goseq(pwf_PSA2vsPSA1,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_PSA2vsPSA0 <- goseq(pwf_PSA2vsPSA0,'hg19','geneSymbol',test.cats="KEGG")

######## PSA1 vs PSA0 ########
# Separo los terminos GO y creo FDRunder y FDRover
GO.BP.PSA1vsPSA0 <- GO.wall_PSA1vsPSA0[GO.wall_PSA1vsPSA0$ontology=="BP",]
GO.BP.PSA1vsPSA0$FDRunder <- p.adjust(GO.BP.PSA1vsPSA0$under_represented_pvalue, n=nrow(GO.BP.PSA1vsPSA0))
GO.BP.PSA1vsPSA0$FDRover <- p.adjust(GO.BP.PSA1vsPSA0$over_represented_pvalue, n=nrow(GO.BP.PSA1vsPSA0))
write.table(GO.BP.PSA1vsPSA0, file = 'GO.biological.process_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)

GO.CC.PSA1vsPSA0 <- GO.wall_PSA1vsPSA0[GO.wall_PSA1vsPSA0$ontology=="CC",]
GO.CC.PSA1vsPSA0$FDRunder <- p.adjust(GO.CC.PSA1vsPSA0$under_represented_pvalue, n=nrow(GO.CC.PSA1vsPSA0))
GO.CC.PSA1vsPSA0$FDRover <- p.adjust(GO.CC.PSA1vsPSA0$over_represented_pvalue, n=nrow(GO.CC.PSA1vsPSA0))
write.table(GO.CC.PSA1vsPSA0, file = 'GO.cellular.component_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)

GO.MF.PSA1vsPSA0 <- GO.wall_PSA1vsPSA0[GO.wall_PSA1vsPSA0$ontology=="MF",]
GO.MF.PSA1vsPSA0$FDRunder <- p.adjust(GO.MF.PSA1vsPSA0$under_represented_pvalue, n=nrow(GO.MF.PSA1vsPSA0))
GO.MF.PSA1vsPSA0$FDRover <- p.adjust(GO.MF.PSA1vsPSA0$over_represented_pvalue, n=nrow(GO.MF.PSA1vsPSA0))
write.table(GO.MF.PSA1vsPSA0, file = 'GO.molecular.function_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.PSA1vsPSA0 <- GO.BP.PSA1vsPSA0[GO.BP.PSA1vsPSA0$FDRover<=0.05,]
enriched.over.GO.CC.PSA1vsPSA0 <- GO.CC.PSA1vsPSA0[GO.CC.PSA1vsPSA0$FDRover<=0.05,]
enriched.over.GO.MF.PSA1vsPSA0 <- GO.MF.PSA1vsPSA0[GO.MF.PSA1vsPSA0$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.PSA1vsPSA0 <- GO.BP.PSA1vsPSA0[GO.BP.PSA1vsPSA0$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.PSA1vsPSA0, file = 'GO.biological.process_under_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.PSA1vsPSA0 <- GO.CC.PSA1vsPSA0[GO.CC.PSA1vsPSA0$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.PSA1vsPSA0, file = 'GO.cellular.component_under_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.PSA1vsPSA0 <- GO.MF.PSA1vsPSA0[GO.MF.PSA1vsPSA0$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.PSA1vsPSA0, file = 'GO.molecular.function_under_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_PSA1vsPSA0$FDRunder <- p.adjust(kegg_DE_PSA1vsPSA0$under_represented_pvalue, n=nrow(kegg_DE_PSA1vsPSA0))
kegg_DE_PSA1vsPSA0$FDRover <- p.adjust(kegg_DE_PSA1vsPSA0$over_represented_pvalue, n=nrow(kegg_DE_PSA1vsPSA0))

kegg_path_PSA1vsPSA0<-mapPathwayToName("hsa")
idx_PSA1vsPSA0<-match(kegg_DE_PSA1vsPSA0$category,kegg_path_PSA1vsPSA0$path)
pathway_name_PSA1vsPSA0<-kegg_path_PSA1vsPSA0[idx_PSA1vsPSA0,2]
data_PSA1vsPSA0<-cbind(kegg_DE_PSA1vsPSA0,pathway_name_PSA1vsPSA0)

# Visualizo el resultado
head(data_PSA1vsPSA0)

write.table(data_PSA1vsPSA0, file = 'KEGG.PATHWAYS_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.PSA1vsPSA0 <- data_PSA1vsPSA0[data_PSA1vsPSA0$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.PSA1vsPSA0 <- data_PSA1vsPSA0[data_PSA1vsPSA0$FDRunder<=0.05,]
write.table(under.KEGG.PSA1vsPSA0, file = 'KEGG.PATHWAYS_under_PSA1vsPSA0.tsv', sep = "\t", row.names = FALSE)

######## PSA2 vs PSA1 ########
# Separo los terminos GO y creo FDRunder y FDRover
GO.BP.PSA2vsPSA1 <- GO.wall_PSA2vsPSA1[GO.wall_PSA2vsPSA1$ontology=="BP",]
GO.BP.PSA2vsPSA1$FDRunder <- p.adjust(GO.BP.PSA2vsPSA1$under_represented_pvalue, n=nrow(GO.BP.PSA2vsPSA1))
GO.BP.PSA2vsPSA1$FDRover <- p.adjust(GO.BP.PSA2vsPSA1$over_represented_pvalue, n=nrow(GO.BP.PSA2vsPSA1))
write.table(GO.BP.PSA2vsPSA1, file = 'GO.biological.process_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)

GO.CC.PSA2vsPSA1 <- GO.wall_PSA2vsPSA1[GO.wall_PSA2vsPSA1$ontology=="CC",]
GO.CC.PSA2vsPSA1$FDRunder <- p.adjust(GO.CC.PSA2vsPSA1$under_represented_pvalue, n=nrow(GO.CC.PSA2vsPSA1))
GO.CC.PSA2vsPSA1$FDRover <- p.adjust(GO.CC.PSA2vsPSA1$over_represented_pvalue, n=nrow(GO.CC.PSA2vsPSA1))
write.table(GO.CC.PSA2vsPSA1, file = 'GO.cellular.component_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)

GO.MF.PSA2vsPSA1 <- GO.wall_PSA2vsPSA1[GO.wall_PSA2vsPSA1$ontology=="MF",]
GO.MF.PSA2vsPSA1$FDRunder <- p.adjust(GO.MF.PSA2vsPSA1$under_represented_pvalue, n=nrow(GO.MF.PSA2vsPSA1))
GO.MF.PSA2vsPSA1$FDRover <- p.adjust(GO.MF.PSA2vsPSA1$over_represented_pvalue, n=nrow(GO.MF.PSA2vsPSA1))
write.table(GO.MF.PSA2vsPSA1, file = 'GO.molecular.function_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.PSA2vsPSA1 <- GO.BP.PSA2vsPSA1[GO.BP.PSA2vsPSA1$FDRover<=0.05,]
enriched.over.GO.CC.PSA2vsPSA1 <- GO.CC.PSA2vsPSA1[GO.CC.PSA2vsPSA1$FDRover<=0.05,]
enriched.over.GO.MF.PSA2vsPSA1 <- GO.MF.PSA2vsPSA1[GO.MF.PSA2vsPSA1$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.PSA2vsPSA1 <- GO.BP.PSA2vsPSA1[GO.BP.PSA2vsPSA1$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.PSA2vsPSA1, file = 'GO.biological.process_under_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.PSA2vsPSA1 <- GO.CC.PSA2vsPSA1[GO.CC.PSA2vsPSA1$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.PSA2vsPSA1, file = 'GO.cellular.component_under_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.PSA2vsPSA1 <- GO.MF.PSA2vsPSA1[GO.MF.PSA2vsPSA1$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.PSA2vsPSA1, file = 'GO.molecular.function_under_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_PSA2vsPSA1$FDRunder <- p.adjust(kegg_DE_PSA2vsPSA1$under_represented_pvalue, n=nrow(kegg_DE_PSA2vsPSA1))
kegg_DE_PSA2vsPSA1$FDRover <- p.adjust(kegg_DE_PSA2vsPSA1$over_represented_pvalue, n=nrow(kegg_DE_PSA2vsPSA1))

kegg_path_PSA2vsPSA1<-mapPathwayToName("hsa")
idx_PSA2vsPSA1<-match(kegg_DE_PSA2vsPSA1$category,kegg_path_PSA2vsPSA1$path)
pathway_name_PSA2vsPSA1<-kegg_path_PSA2vsPSA1[idx_PSA2vsPSA1,2]
data_PSA2vsPSA1<-cbind(kegg_DE_PSA2vsPSA1,pathway_name_PSA2vsPSA1)

# Visualizo el resultado
head(data_PSA2vsPSA1)

write.table(data_PSA2vsPSA1, file = 'KEGG.PATHWAYS_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.PSA2vsPSA1 <- data_PSA2vsPSA1[data_PSA2vsPSA1$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.PSA2vsPSA1 <- data_PSA2vsPSA1[data_PSA2vsPSA1$FDRunder<=0.05,]
write.table(under.KEGG.PSA2vsPSA1, file = 'KEGG.PATHWAYS_under_PSA2vsPSA1.tsv', sep = "\t", row.names = FALSE)


######## PSA2 vs PSA0 ########
# Separo los terminos GO y creo FDRdown y FDRover
GO.BP.PSA2vsPSA0 <- GO.wall_PSA2vsPSA0[GO.wall_PSA2vsPSA0$ontology=="BP",]
GO.BP.PSA2vsPSA0$FDRunder <- p.adjust(GO.BP.PSA2vsPSA0$under_represented_pvalue, n=nrow(GO.BP.PSA2vsPSA0))
GO.BP.PSA2vsPSA0$FDRover <- p.adjust(GO.BP.PSA2vsPSA0$over_represented_pvalue, n=nrow(GO.BP.PSA2vsPSA0))
write.table(GO.BP.PSA2vsPSA0, file = 'GO.biological.process_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)

GO.CC.PSA2vsPSA0 <- GO.wall_PSA2vsPSA0[GO.wall_PSA2vsPSA0$ontology=="CC",]
GO.CC.PSA2vsPSA0$FDRunder <- p.adjust(GO.CC.PSA2vsPSA0$under_represented_pvalue, n=nrow(GO.CC.PSA2vsPSA0))
GO.CC.PSA2vsPSA0$FDRover <- p.adjust(GO.CC.PSA2vsPSA0$over_represented_pvalue, n=nrow(GO.CC.PSA2vsPSA0))
write.table(GO.CC.PSA2vsPSA0, file = 'GO.cellular.component_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)

GO.MF.PSA2vsPSA0 <- GO.wall_PSA2vsPSA0[GO.wall_PSA2vsPSA0$ontology=="MF",]
GO.MF.PSA2vsPSA0$FDRunder <- p.adjust(GO.MF.PSA2vsPSA0$under_represented_pvalue, n=nrow(GO.MF.PSA2vsPSA0))
GO.MF.PSA2vsPSA0$FDRover <- p.adjust(GO.MF.PSA2vsPSA0$over_represented_pvalue, n=nrow(GO.MF.PSA2vsPSA0))
write.table(GO.MF.PSA2vsPSA0, file = 'GO.molecular.function_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.PSA2vsPSA0 <- GO.BP.PSA2vsPSA0[GO.BP.PSA2vsPSA0$FDRover<=0.05,]
enriched.over.GO.CC.PSA2vsPSA0 <- GO.CC.PSA2vsPSA0[GO.CC.PSA2vsPSA0$FDRover<=0.05,]
enriched.over.GO.MF.PSA2vsPSA0 <- GO.MF.PSA2vsPSA0[GO.MF.PSA2vsPSA0$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.PSA2vsPSA0 <- GO.BP.PSA2vsPSA0[GO.BP.PSA2vsPSA0$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.PSA2vsPSA0, file = 'GO.biological.process_under_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.PSA2vsPSA0 <- GO.CC.PSA2vsPSA0[GO.CC.PSA2vsPSA0$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.PSA2vsPSA0, file = 'GO.cellular.component_under_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.PSA2vsPSA0 <- GO.MF.PSA2vsPSA0[GO.MF.PSA2vsPSA0$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.PSA2vsPSA0, file = 'GO.molecular.function_under_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_PSA2vsPSA0$FDRunder <- p.adjust(kegg_DE_PSA2vsPSA0$under_represented_pvalue, n=nrow(kegg_DE_PSA2vsPSA0))
kegg_DE_PSA2vsPSA0$FDRover <- p.adjust(kegg_DE_PSA2vsPSA0$over_represented_pvalue, n=nrow(kegg_DE_PSA2vsPSA0))

kegg_path_PSA2vsPSA0<-mapPathwayToName("hsa")
idx_PSA2vsPSA0<-match(kegg_DE_PSA2vsPSA0$category,kegg_path_PSA2vsPSA0$path)
pathway_name_PSA2vsPSA0<-kegg_path_PSA2vsPSA0[idx_PSA2vsPSA0,2]
data_PSA2vsPSA0<-cbind(kegg_DE_PSA2vsPSA0,pathway_name_PSA2vsPSA0)

# Visualizo el resultado
head(data_PSA2vsPSA0)

write.table(data_PSA2vsPSA0, file = 'KEGG.PATHWAYS_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.PSA2vsPSA0 <- data_PSA2vsPSA0[data_PSA2vsPSA0$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.PSA2vsPSA0 <- data_PSA2vsPSA0[data_PSA2vsPSA0$FDRunder<=0.05,]
write.table(under.KEGG.PSA2vsPSA0, file = 'KEGG.PATHWAYS_under_PSA2vsPSA0.tsv', sep = "\t", row.names = FALSE)