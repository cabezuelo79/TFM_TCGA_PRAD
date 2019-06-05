##############################################################
####### Analisis de expresion genica mediante RNA-Seq#########        
##############       GLEASON   ###############################
###### Analisis de expresion diferencial con el paquete edgeR#
##############################################################

# Establecimiento de directorio de trabajo, donde se encuentra el archivo "GENERAL.txt"
setwd('C:/Users/JoseMaria/OneDrive - MERIDIEM SEEDS S.L/Documentos/Master/TFM/RNASeq_PRAD/GitHub/TFM_TCGA_PRAD')

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

# Cargo el vector 'Gleason' de forma manual
Gleason <- factor(c("G1","G1","G2","G1","G2","G2","G2","G2","G1","G1","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G2","G2","G2","G2","G2","G0","G2","G0","G1","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G0","G2","G0","G2","G2","G2","G0","G2","G0","G2","G0","G2","G2","G0","G1","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G0","G2","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G1","G2","G2","G0","G2","G2","G2","G2","G0","G1","G2","G1","G0","G2","G0","G2","G2","G2","G0","G2","G2","G0","G2","G0","G2","G2","G0","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G0","G2","G2","G2","G0","G2","G2","G0","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G1","G1","G2","G2","G2","G2","G1","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G0","G2","G0","G2","G1","G2","G2","G0","G2","G2","G2","G0","G2","G2","G2","G1","G2","G2","G2","G0","G1","G0","G1","G0","G2","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G0","G2","G1","G2","G2","G2","G1","G2","G1","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G1","G2","G2","G1","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2"
))

#Obtengo el tamaño de cada gen para normalizar teniendo en ceunta del tamaño de los genes
gene.length_G <-read.table("his-Size.tab",header=T)
idx_G <-match(rawdata$SYMBOL,gene.length_G$Gene)
results_counts_G <-gene.length_G[idx_G,]
results_counts_G[is.na(results_counts_G$Length),"Length"]<-0
nrow(results_counts_G)

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes, el resto las muestras
# como grupo selecciono tipo de tejido 'Gleason'
yG <- DGEList(counts=rawdata[,3:551], genes=results_counts_G, group=Gleason)

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
top_G1vsG0_filter <- as.vector(na.omit(top_G1vsG0$table[top_G1vsG0$table$FDR < 0.05,"Gene"]))
top_G2vsG1_filter <- as.vector(na.omit(top_G2vsG1$table[top_G2vsG1$table$FDR < 0.05,"Gene"]))
top_G2vsG0_filter <- as.vector(na.omit(top_G2vsG0$table[top_G2vsG0$table$FDR < 0.05,"Gene"]))

# Genero vectores 'universal' que contienen tanto los genes DE como los no DE
universeG1vsG0 <- as.vector(unique(na.omit(top_G1vsG0$table$Gene)))
universeG2vsG1 <- as.vector(unique(na.omit(top_G2vsG1$table$Gene)))
universeG2vsG0 <- as.vector(unique(na.omit(top_G2vsG0$table$Gene)))

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


# Realizo los analisis GO mediante la aproximación de Wallenius
GO.wall_G1vsG0=goseq(pwf_G1vsG0,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
GO.wall_G2vsG1=goseq(pwf_G2vsG1,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
GO.wall_G2vsG0=goseq(pwf_G2vsG0,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
kegg_DE_G1vsG0 <- goseq(pwf_G1vsG0,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_G2vsG1 <- goseq(pwf_G2vsG1,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_G2vsG0 <- goseq(pwf_G2vsG0,'hg19','geneSymbol',test.cats="KEGG")

######## G1 vs G0 ########
# Separo los terminos GO y creo FDRunder y FDRover
GO.BP.G1vsG0 <- GO.wall_G1vsG0[GO.wall_G1vsG0$ontology=="BP",]
GO.BP.G1vsG0$FDRunder <- p.adjust(GO.BP.G1vsG0$under_represented_pvalue, n=nrow(GO.BP.G1vsG0))
GO.BP.G1vsG0$FDRover <- p.adjust(GO.BP.G1vsG0$over_represented_pvalue, n=nrow(GO.BP.G1vsG0))
write.table(GO.BP.G1vsG0, file = 'GO.biological.process_G1vsG0.tsv', sep = "\t", row.names = FALSE)

GO.CC.G1vsG0 <- GO.wall_G1vsG0[GO.wall_G1vsG0$ontology=="CC",]
GO.CC.G1vsG0$FDRunder <- p.adjust(GO.CC.G1vsG0$under_represented_pvalue, n=nrow(GO.CC.G1vsG0))
GO.CC.G1vsG0$FDRover <- p.adjust(GO.CC.G1vsG0$over_represented_pvalue, n=nrow(GO.CC.G1vsG0))
write.table(GO.CC.G1vsG0, file = 'GO.cellular.component_G1vsG0.tsv', sep = "\t", row.names = FALSE)

GO.MF.G1vsG0 <- GO.wall_G1vsG0[GO.wall_G1vsG0$ontology=="MF",]
GO.MF.G1vsG0$FDRunder <- p.adjust(GO.MF.G1vsG0$under_represented_pvalue, n=nrow(GO.MF.G1vsG0))
GO.MF.G1vsG0$FDRover <- p.adjust(GO.MF.G1vsG0$over_represented_pvalue, n=nrow(GO.MF.G1vsG0))
write.table(GO.MF.G1vsG0, file = 'GO.molecular.function_G1vsG0.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.G1vsG0 <- GO.BP.G1vsG0[GO.BP.G1vsG0$FDRover<=0.05,]
enriched.over.GO.CC.G1vsG0 <- GO.CC.G1vsG0[GO.CC.G1vsG0$FDRover<=0.05,]
enriched.over.GO.MF.G1vsG0 <- GO.MF.G1vsG0[GO.MF.G1vsG0$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.G1vsG0 <- GO.BP.G1vsG0[GO.BP.G1vsG0$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.G1vsG0, file = 'GO.biological.process_under_G1vsG0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.G1vsG0 <- GO.CC.G1vsG0[GO.CC.G1vsG0$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.G1vsG0, file = 'GO.cellular.component_under_G1vsG0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.G1vsG0 <- GO.MF.G1vsG0[GO.MF.G1vsG0$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.G1vsG0, file = 'GO.molecular.function_under_G1vsG0.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_G1vsG0$FDRunder <- p.adjust(kegg_DE_G1vsG0$under_represented_pvalue, n=nrow(kegg_DE_G1vsG0))
kegg_DE_G1vsG0$FDRover <- p.adjust(kegg_DE_G1vsG0$over_represented_pvalue, n=nrow(kegg_DE_G1vsG0))

kegg_path_G1vsG0<-mapPathwayToName("hsa")
idx_G1vsG0<-match(kegg_DE_G1vsG0$category,kegg_path_G1vsG0$path)
pathway_name_G1vsG0<-kegg_path_G1vsG0[idx_G1vsG0,2]
data_G1vsG0<-cbind(kegg_DE_G1vsG0,pathway_name_G1vsG0)

# Visualizo el resultado
head(data_G1vsG0)

write.table(data_G1vsG0, file = 'KEGG.PATHWAYS_G1vsG0.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.G1vsG0 <- data_G1vsG0[data_G1vsG0$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.G1vsG0 <- data_G1vsG0[data_G1vsG0$FDRunder<=0.05,]
write.table(under.KEGG.G1vsG0, file = 'KEGG.PATHWAYS_under_G1vsG0.tsv', sep = "\t", row.names = FALSE)

######## G2 vs G1 ########
# Separo los terminos GO y creo FDRunder y FDRover
GO.BP.G2vsG1 <- GO.wall_G2vsG1[GO.wall_G2vsG1$ontology=="BP",]
GO.BP.G2vsG1$FDRunder <- p.adjust(GO.BP.G2vsG1$under_represented_pvalue, n=nrow(GO.BP.G2vsG1))
GO.BP.G2vsG1$FDRover <- p.adjust(GO.BP.G2vsG1$over_represented_pvalue, n=nrow(GO.BP.G2vsG1))
write.table(GO.BP.G2vsG1, file = 'GO.biological.process_G2vsG1.tsv', sep = "\t", row.names = FALSE)

GO.CC.G2vsG1 <- GO.wall_G2vsG1[GO.wall_G2vsG1$ontology=="CC",]
GO.CC.G2vsG1$FDRunder <- p.adjust(GO.CC.G2vsG1$under_represented_pvalue, n=nrow(GO.CC.G2vsG1))
GO.CC.G2vsG1$FDRover <- p.adjust(GO.CC.G2vsG1$over_represented_pvalue, n=nrow(GO.CC.G2vsG1))
write.table(GO.CC.G2vsG1, file = 'GO.cellular.component_G2vsG1.tsv', sep = "\t", row.names = FALSE)

GO.MF.G2vsG1 <- GO.wall_G2vsG1[GO.wall_G2vsG1$ontology=="MF",]
GO.MF.G2vsG1$FDRunder <- p.adjust(GO.MF.G2vsG1$under_represented_pvalue, n=nrow(GO.MF.G2vsG1))
GO.MF.G2vsG1$FDRover <- p.adjust(GO.MF.G2vsG1$over_represented_pvalue, n=nrow(GO.MF.G2vsG1))
write.table(GO.MF.G2vsG1, file = 'GO.molecular.function_G2vsG1.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.G2vsG1 <- GO.BP.G2vsG1[GO.BP.G2vsG1$FDRover<=0.05,]
enriched.over.GO.CC.G2vsG1 <- GO.CC.G2vsG1[GO.CC.G2vsG1$FDRover<=0.05,]
enriched.over.GO.MF.G2vsG1 <- GO.MF.G2vsG1[GO.MF.G2vsG1$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.G2vsG1 <- GO.BP.G2vsG1[GO.BP.G2vsG1$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.G2vsG1, file = 'GO.biological.process_under_G2vsG1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.G2vsG1 <- GO.CC.G2vsG1[GO.CC.G2vsG1$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.G2vsG1, file = 'GO.cellular.component_under_G2vsG1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.G2vsG1 <- GO.MF.G2vsG1[GO.MF.G2vsG1$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.G2vsG1, file = 'GO.molecular.function_under_G2vsG1.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_G2vsG1$FDRunder <- p.adjust(kegg_DE_G2vsG1$under_represented_pvalue, n=nrow(kegg_DE_G2vsG1))
kegg_DE_G2vsG1$FDRover <- p.adjust(kegg_DE_G2vsG1$over_represented_pvalue, n=nrow(kegg_DE_G2vsG1))

kegg_path_G2vsG1<-mapPathwayToName("hsa")
idx_G2vsG1<-match(kegg_DE_G2vsG1$category,kegg_path_G2vsG1$path)
pathway_name_G2vsG1<-kegg_path_G2vsG1[idx_G2vsG1,2]
data_G2vsG1<-cbind(kegg_DE_G2vsG1,pathway_name_G2vsG1)

# Visualizo el resultado
head(data_G2vsG1)

write.table(data_G2vsG1, file = 'KEGG.PATHWAYS_G2vsG1.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.G2vsG1 <- data_G2vsG1[data_G2vsG1$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.G2vsG1 <- data_G2vsG1[data_G2vsG1$FDRunder<=0.05,]
write.table(under.KEGG.G2vsG1, file = 'KEGG.PATHWAYS_under_G2vsG1.tsv', sep = "\t", row.names = FALSE)


######## G2 vs G0 ########
# Separo los terminos GO y creo FDRdown y FDRover
GO.BP.G2vsG0 <- GO.wall_G2vsG0[GO.wall_G2vsG0$ontology=="BP",]
GO.BP.G2vsG0$FDRunder <- p.adjust(GO.BP.G2vsG0$under_represented_pvalue, n=nrow(GO.BP.G2vsG0))
GO.BP.G2vsG0$FDRover <- p.adjust(GO.BP.G2vsG0$over_represented_pvalue, n=nrow(GO.BP.G2vsG0))
write.table(GO.BP.G2vsG0, file = 'GO.biological.process_G2vsG0.tsv', sep = "\t", row.names = FALSE)

GO.CC.G2vsG0 <- GO.wall_G2vsG0[GO.wall_G2vsG0$ontology=="CC",]
GO.CC.G2vsG0$FDRunder <- p.adjust(GO.CC.G2vsG0$under_represented_pvalue, n=nrow(GO.CC.G2vsG0))
GO.CC.G2vsG0$FDRover <- p.adjust(GO.CC.G2vsG0$over_represented_pvalue, n=nrow(GO.CC.G2vsG0))
write.table(GO.CC.G2vsG0, file = 'GO.cellular.component_G2vsG0.tsv', sep = "\t", row.names = FALSE)

GO.MF.G2vsG0 <- GO.wall_G2vsG0[GO.wall_G2vsG0$ontology=="MF",]
GO.MF.G2vsG0$FDRunder <- p.adjust(GO.MF.G2vsG0$under_represented_pvalue, n=nrow(GO.MF.G2vsG0))
GO.MF.G2vsG0$FDRover <- p.adjust(GO.MF.G2vsG0$over_represented_pvalue, n=nrow(GO.MF.G2vsG0))
write.table(GO.MF.G2vsG0, file = 'GO.molecular.function_G2vsG0.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.G2vsG0 <- GO.BP.G2vsG0[GO.BP.G2vsG0$FDRover<=0.05,]
enriched.over.GO.CC.G2vsG0 <- GO.CC.G2vsG0[GO.CC.G2vsG0$FDRover<=0.05,]
enriched.over.GO.MF.G2vsG0 <- GO.MF.G2vsG0[GO.MF.G2vsG0$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.G2vsG0 <- GO.BP.G2vsG0[GO.BP.G2vsG0$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.G2vsG0, file = 'GO.biological.process_under_G2vsG0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.G2vsG0 <- GO.CC.G2vsG0[GO.CC.G2vsG0$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.G2vsG0, file = 'GO.cellular.component_under_G2vsG0.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.G2vsG0 <- GO.MF.G2vsG0[GO.MF.G2vsG0$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.G2vsG0, file = 'GO.molecular.function_under_G2vsG0.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_G2vsG0$FDRunder <- p.adjust(kegg_DE_G2vsG0$under_represented_pvalue, n=nrow(kegg_DE_G2vsG0))
kegg_DE_G2vsG0$FDRover <- p.adjust(kegg_DE_G2vsG0$over_represented_pvalue, n=nrow(kegg_DE_G2vsG0))

kegg_path_G2vsG0<-mapPathwayToName("hsa")
idx_G2vsG0<-match(kegg_DE_G2vsG0$category,kegg_path_G2vsG0$path)
pathway_name_G2vsG0<-kegg_path_G2vsG0[idx_G2vsG0,2]
data_G2vsG0<-cbind(kegg_DE_G2vsG0,pathway_name_G2vsG0)

# Visualizo el resultado
head(data_G2vsG0)

write.table(data_G2vsG0, file = 'KEGG.PATHWAYS_G2vsG0.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.G2vsG0 <- data_G2vsG0[data_G2vsG0$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.G2vsG0 <- data_G2vsG0[data_G2vsG0$FDRunder<=0.05,]
write.table(under.KEGG.G2vsG0, file = 'KEGG.PATHWAYS_under_G2vsG0.tsv', sep = "\t", row.names = FALSE)