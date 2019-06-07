##############################################################
####### Analisis de expresion genica mediante RNA-Seq#########        
##############       DRUGS     ###############################
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
rawdata_drugs <- read.delim("DRUGS.txt", check.names=FALSE, stringsAsFactors=FALSE)

#incluyo el nombre de los genes como rowname
rownames(rawdata_drugs)<-rawdata_drugs$SYMBOL

# Cargo el vector 'Drugs' de forma manual
Drugs <- factor(c("DR1","DR1","DR1","DR2","DR1","DR1","DR1","DR1","DR2","DR3","DR1","DR1","DR1","DR1","DR1","DR1","DR1","DR2","DR1","DR3","DR1","DR1","DR1","DR1","DR3","DR1","DR2","DR1","DR2","DR2","DR1","DR2","DR1","DR1","DR1","DR1","DR2","DR1","DR3","DR2","DR1","DR1","DR1","DR1","DR1","DR1","DR2","DR1","DR1","DR3","DR1","DR2","DR2","DR2","DR2","DR2","DR2","DR2","DR2","DR1","DR1","DR1","DR2","DR2","DR2","DR2","DR2","DR2","DR1","DR2","DR1","DR1","DR1"
))

#Obtengo el tamaño de cada gen para normalizar teniendo en ceunta del tamaño de los genes
gene.length_D <-read.table("his-Size.tab",header=T)
idx_D <-match(rawdata$SYMBOL,gene.length_D$Gene)
results_counts_D <-gene.length_D[idx_D,]
results_counts_D[is.na(results_counts_D$Length),"Length"]<-0
nrow(results_counts_D)

# Construyo un objeto DGEList con la información del archivo, las dos primeras columnas son el simbolo e identidad de los genes, el resto las muestras
# como grupo selecciono tipo de tratamiento 'Drugs'
yD <- DGEList(counts=rawdata_drugs[,3:75], genes=results_counts_D, group=Drugs)

# Llevo a cabo la normalizacion de las muestras
yD <- calcNormFactors(yD)

#Diseño la matriz para glm
designD <- model.matrix(~0+Drugs, data=yD$samples)
colnames(designD) <- levels(yD$samples$group)

# visualizo la matriz
designD

# Calculo las dispersiones
yD <- estimateGLMCommonDisp(yD, designD)
yD <- estimateGLMTrendedDisp(yD, designD)
yD <- estimateGLMTagwiseDisp(yD, designD)

# visualizo la dispersion en un grafico BCV
plotBCV(yD)

# Quasi-likelihood test
fitD <- glmQLFit(yD, designD)

#CONTRASTES
my.contrasts_D <- makeContrasts(DR2vsDR1=DR2-DR1, DR3vsDR2=DR3-DR2, DR3vsDR1=DR3-DR1, levels=designD)
qlf.DR2vsDR1 <- glmQLFTest(fitD, contrast=my.contrasts_D[,"DR2vsDR1"])
top_DR2vsDR1 <- topTags(qlf.DR2vsDR1, n=Inf)
qlf.DR3vsDR2 <- glmQLFTest(fitD, contrast=my.contrasts_D[,"DR3vsDR2"])
top_DR3vsDR2 <- topTags(qlf.DR3vsDR2, n=Inf)
qlf.DR3vsDR1 <- glmQLFTest(fitD, contrast=my.contrasts_D[,"DR3vsDR1"])
top_DR3vsDR1 <- topTags(qlf.DR3vsDR1, n=Inf)

# visualizo el resumen de genes sobre-expresados, genes sub-expresados y genes no significativos
summary(decideTests(qlf.DR2vsDR1))
summary(decideTests(qlf.DR3vsDR2))
summary(decideTests(qlf.DR3vsDR1))

##############################################################################################
############## ANALISIS DE ENRIQUECIMIENTO FUNCIONAL GOSEQ      ##############################
##############################################################################################

# Cargo los paquetes necesarios para el analisis
library(goseq)
library(GO.db)

# Genero vectores con los genes diferencialmente expresados (FDR < 0.05)
top_DR2vsDR1_filter <- as.vector(na.omit(top_DR2vsDR1$table[top_DR2vsDR1$table$FDR < 0.05,"Gene"]))
top_DR3vsDR2_filter <- as.vector(na.omit(top_DR3vsDR2$table[top_DR3vsDR2$table$FDR < 0.05,"Gene"]))
top_DR3vsDR1_filter <- as.vector(na.omit(top_DR3vsDR1$table[top_DR3vsDR1$table$FDR < 0.05,"Gene"]))

# Genero vectores 'universal' que contienen tanto los genes DE como los no DE
universeDR2vsDR1 <- as.vector(unique(na.omit(top_DR2vsDR1$table$Gene)))
universeDR3vsDR2 <- as.vector(unique(na.omit(top_DR3vsDR2$table$Gene)))
universeDR3vsDR1 <- as.vector(unique(na.omit(top_DR3vsDR1$table$Gene)))

# Construyo vectores apropiados para el analisis con 'goseq'
mytableDR2vsDR1 <- as.integer(unique(universeDR2vsDR1)%in%top_DR2vsDR1_filter)
names(mytableDR2vsDR1) <- unique(universeDR2vsDR1)
mytableDR3vsDR2 <- as.integer(unique(universeDR3vsDR2)%in%top_DR3vsDR2_filter)
names(mytableDR3vsDR2) <- unique(universeDR3vsDR2)
mytableDR3vsDR1 <- as.integer(unique(universeDR3vsDR1)%in%top_DR3vsDR1_filter)
names(mytableDR3vsDR1) <- unique(universeDR3vsDR1)

# Veo la cantidad de genes DE (1) y no DE (0) y compruebo que este correcto
table(mytableDR2vsDR1)
head(mytableDR2vsDR1)
table(mytableDR3vsDR2)
head(mytableDR3vsDR2)
table(mytableDR3vsDR1)
head(mytableDR3vsDR1)


# Probability Weighting Function (pwf) utilizando 'geneSymbol' y genoma 'hg19'
pwf_DR2vsDR1=nullp(mytableDR2vsDR1,"hg19","geneSymbol")
pwf_DR3vsDR2=nullp(mytableDR3vsDR2,"hg19","geneSymbol")
pwf_DR3vsDR1=nullp(mytableDR3vsDR1,"hg19","geneSymbol")

# Compruebo el resultado
head(pwf_DR2vsDR1)
head(pwf_DR3vsDR2)
head(pwf_DR3vsDR1)


# Realizo los analisis GO mediante la aproximación de Wallenius
GO.wall_DR2vsDR1=goseq(pwf_DR2vsDR1,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
GO.wall_DR3vsDR2=goseq(pwf_DR3vsDR2,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
GO.wall_DR3vsDR1=goseq(pwf_DR3vsDR1,"hg19","geneSymbol",test.cats = c("GO:BP","GO:CC","GO:MF"))
kegg_DE_DR2vsDR1 <- goseq(pwf_DR2vsDR1,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_DR3vsDR2 <- goseq(pwf_DR3vsDR2,'hg19','geneSymbol',test.cats="KEGG")
kegg_DE_DR3vsDR1 <- goseq(pwf_DR3vsDR1,'hg19','geneSymbol',test.cats="KEGG")

######## DR2 vs DR1 ########
# Separo los terminos GO y creo FDRunder y FDRover
GO.BP.DR2vsDR1 <- GO.wall_DR2vsDR1[GO.wall_DR2vsDR1$ontology=="BP",]
GO.BP.DR2vsDR1$FDRunder <- p.adjust(GO.BP.DR2vsDR1$under_represented_pvalue, n=nrow(GO.BP.DR2vsDR1))
GO.BP.DR2vsDR1$FDRover <- p.adjust(GO.BP.DR2vsDR1$over_represented_pvalue, n=nrow(GO.BP.DR2vsDR1))
write.table(GO.BP.DR2vsDR1, file = 'GO.biological.process_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)

GO.CC.DR2vsDR1 <- GO.wall_DR2vsDR1[GO.wall_DR2vsDR1$ontology=="CC",]
GO.CC.DR2vsDR1$FDRunder <- p.adjust(GO.CC.DR2vsDR1$under_represented_pvalue, n=nrow(GO.CC.DR2vsDR1))
GO.CC.DR2vsDR1$FDRover <- p.adjust(GO.CC.DR2vsDR1$over_represented_pvalue, n=nrow(GO.CC.DR2vsDR1))
write.table(GO.CC.DR2vsDR1, file = 'GO.cellular.component_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)

GO.MF.DR2vsDR1 <- GO.wall_DR2vsDR1[GO.wall_DR2vsDR1$ontology=="MF",]
GO.MF.DR2vsDR1$FDRunder <- p.adjust(GO.MF.DR2vsDR1$under_represented_pvalue, n=nrow(GO.MF.DR2vsDR1))
GO.MF.DR2vsDR1$FDRover <- p.adjust(GO.MF.DR2vsDR1$over_represented_pvalue, n=nrow(GO.MF.DR2vsDR1))
write.table(GO.MF.DR2vsDR1, file = 'GO.molecular.function_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.DR2vsDR1 <- GO.BP.DR2vsDR1[GO.BP.DR2vsDR1$FDRover<=0.05,]
enriched.over.GO.CC.DR2vsDR1 <- GO.CC.DR2vsDR1[GO.CC.DR2vsDR1$FDRover<=0.05,]
enriched.over.GO.MF.DR2vsDR1 <- GO.MF.DR2vsDR1[GO.MF.DR2vsDR1$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.DR2vsDR1 <- GO.BP.DR2vsDR1[GO.BP.DR2vsDR1$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.DR2vsDR1, file = 'GO.biological.process_under_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.DR2vsDR1 <- GO.CC.DR2vsDR1[GO.CC.DR2vsDR1$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.DR2vsDR1, file = 'GO.cellular.component_under_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.DR2vsDR1 <- GO.MF.DR2vsDR1[GO.MF.DR2vsDR1$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.DR2vsDR1, file = 'GO.molecular.function_under_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_DR2vsDR1$FDRunder <- p.adjust(kegg_DE_DR2vsDR1$under_represented_pvalue, n=nrow(kegg_DE_DR2vsDR1))
kegg_DE_DR2vsDR1$FDRover <- p.adjust(kegg_DE_DR2vsDR1$over_represented_pvalue, n=nrow(kegg_DE_DR2vsDR1))

kegg_path_DR2vsDR1<-mapPathwayToName("hsa")
idx_DR2vsDR1<-match(kegg_DE_DR2vsDR1$category,kegg_path_DR2vsDR1$path)
pathway_name_DR2vsDR1<-kegg_path_DR2vsDR1[idx_DR2vsDR1,2]
data_DR2vsDR1<-cbind(kegg_DE_DR2vsDR1,pathway_name_DR2vsDR1)

# Visualizo el resultado
head(data_DR2vsDR1)

write.table(data_DR2vsDR1, file = 'KEGG.PATHWAYS_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.DR2vsDR1 <- data_DR2vsDR1[data_DR2vsDR1$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.DR2vsDR1 <- data_DR2vsDR1[data_DR2vsDR1$FDRunder<=0.05,]
write.table(under.KEGG.DR2vsDR1, file = 'KEGG.PATHWAYS_under_DR2vsDR1.tsv', sep = "\t", row.names = FALSE)

######## DR3 vs DR2 ########
# Separo los terminos GO y creo FDRunder y FDRover
GO.BP.DR3vsDR2 <- GO.wall_DR3vsDR2[GO.wall_DR3vsDR2$ontology=="BP",]
GO.BP.DR3vsDR2$FDRunder <- p.adjust(GO.BP.DR3vsDR2$under_represented_pvalue, n=nrow(GO.BP.DR3vsDR2))
GO.BP.DR3vsDR2$FDRover <- p.adjust(GO.BP.DR3vsDR2$over_represented_pvalue, n=nrow(GO.BP.DR3vsDR2))
write.table(GO.BP.DR3vsDR2, file = 'GO.biological.process_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)

GO.CC.DR3vsDR2 <- GO.wall_DR3vsDR2[GO.wall_DR3vsDR2$ontology=="CC",]
GO.CC.DR3vsDR2$FDRunder <- p.adjust(GO.CC.DR3vsDR2$under_represented_pvalue, n=nrow(GO.CC.DR3vsDR2))
GO.CC.DR3vsDR2$FDRover <- p.adjust(GO.CC.DR3vsDR2$over_represented_pvalue, n=nrow(GO.CC.DR3vsDR2))
write.table(GO.CC.DR3vsDR2, file = 'GO.cellular.component_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)

GO.MF.DR3vsDR2 <- GO.wall_DR3vsDR2[GO.wall_DR3vsDR2$ontology=="MF",]
GO.MF.DR3vsDR2$FDRunder <- p.adjust(GO.MF.DR3vsDR2$under_represented_pvalue, n=nrow(GO.MF.DR3vsDR2))
GO.MF.DR3vsDR2$FDRover <- p.adjust(GO.MF.DR3vsDR2$over_represented_pvalue, n=nrow(GO.MF.DR3vsDR2))
write.table(GO.MF.DR3vsDR2, file = 'GO.molecular.function_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.DR3vsDR2 <- GO.BP.DR3vsDR2[GO.BP.DR3vsDR2$FDRover<=0.05,]
enriched.over.GO.CC.DR3vsDR2 <- GO.CC.DR3vsDR2[GO.CC.DR3vsDR2$FDRover<=0.05,]
enriched.over.GO.MF.DR3vsDR2 <- GO.MF.DR3vsDR2[GO.MF.DR3vsDR2$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.DR3vsDR2 <- GO.BP.DR3vsDR2[GO.BP.DR3vsDR2$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.DR3vsDR2, file = 'GO.biological.process_under_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.DR3vsDR2 <- GO.CC.DR3vsDR2[GO.CC.DR3vsDR2$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.DR3vsDR2, file = 'GO.cellular.component_under_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.DR3vsDR2 <- GO.MF.DR3vsDR2[GO.MF.DR3vsDR2$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.DR3vsDR2, file = 'GO.molecular.function_under_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_DR3vsDR2$FDRunder <- p.adjust(kegg_DE_DR3vsDR2$under_represented_pvalue, n=nrow(kegg_DE_DR3vsDR2))
kegg_DE_DR3vsDR2$FDRover <- p.adjust(kegg_DE_DR3vsDR2$over_represented_pvalue, n=nrow(kegg_DE_DR3vsDR2))

kegg_path_DR3vsDR2<-mapPathwayToName("hsa")
idx_DR3vsDR2<-match(kegg_DE_DR3vsDR2$category,kegg_path_DR3vsDR2$path)
pathway_name_DR3vsDR2<-kegg_path_DR3vsDR2[idx_DR3vsDR2,2]
data_DR3vsDR2<-cbind(kegg_DE_DR3vsDR2,pathway_name_DR3vsDR2)

# Visualizo el resultado
head(data_DR3vsDR2)

write.table(data_DR3vsDR2, file = 'KEGG.PATHWAYS_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.DR3vsDR2 <- data_DR3vsDR2[data_DR3vsDR2$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.DR3vsDR2 <- data_DR3vsDR2[data_DR3vsDR2$FDRunder<=0.05,]
write.table(under.KEGG.DR3vsDR2, file = 'KEGG.PATHWAYS_under_DR3vsDR2.tsv', sep = "\t", row.names = FALSE)


######## DR3 vs DR1 ########
# Separo los terminos GO y creo FDRdown y FDRover
GO.BP.DR3vsDR1 <- GO.wall_DR3vsDR1[GO.wall_DR3vsDR1$ontology=="BP",]
GO.BP.DR3vsDR1$FDRunder <- p.adjust(GO.BP.DR3vsDR1$under_represented_pvalue, n=nrow(GO.BP.DR3vsDR1))
GO.BP.DR3vsDR1$FDRover <- p.adjust(GO.BP.DR3vsDR1$over_represented_pvalue, n=nrow(GO.BP.DR3vsDR1))
write.table(GO.BP.DR3vsDR1, file = 'GO.biological.process_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)

GO.CC.DR3vsDR1 <- GO.wall_DR3vsDR1[GO.wall_DR3vsDR1$ontology=="CC",]
GO.CC.DR3vsDR1$FDRunder <- p.adjust(GO.CC.DR3vsDR1$under_represented_pvalue, n=nrow(GO.CC.DR3vsDR1))
GO.CC.DR3vsDR1$FDRover <- p.adjust(GO.CC.DR3vsDR1$over_represented_pvalue, n=nrow(GO.CC.DR3vsDR1))
write.table(GO.CC.DR3vsDR1, file = 'GO.cellular.component_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)

GO.MF.DR3vsDR1 <- GO.wall_DR3vsDR1[GO.wall_DR3vsDR1$ontology=="MF",]
GO.MF.DR3vsDR1$FDRunder <- p.adjust(GO.MF.DR3vsDR1$under_represented_pvalue, n=nrow(GO.MF.DR3vsDR1))
GO.MF.DR3vsDR1$FDRover <- p.adjust(GO.MF.DR3vsDR1$over_represented_pvalue, n=nrow(GO.MF.DR3vsDR1))
write.table(GO.MF.DR3vsDR1, file = 'GO.molecular.function_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)

# Guardo los FDR <= 0.05
enriched.over.GO.BP.DR3vsDR1 <- GO.BP.DR3vsDR1[GO.BP.DR3vsDR1$FDRover<=0.05,]
enriched.over.GO.CC.DR3vsDR1 <- GO.CC.DR3vsDR1[GO.CC.DR3vsDR1$FDRover<=0.05,]
enriched.over.GO.MF.DR3vsDR1 <- GO.MF.DR3vsDR1[GO.MF.DR3vsDR1$FDRover<=0.05,]

# Genes under-expressed
enriched.under.GO.BP.DR3vsDR1 <- GO.BP.DR3vsDR1[GO.BP.DR3vsDR1$FDRunder<=0.05,]
write.table(enriched.under.GO.BP.DR3vsDR1, file = 'GO.biological.process_under_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.CC.DR3vsDR1 <- GO.CC.DR3vsDR1[GO.CC.DR3vsDR1$FDRunder<=0.05,]
write.table(enriched.under.GO.CC.DR3vsDR1, file = 'GO.cellular.component_under_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)
enriched.under.GO.MF.DR3vsDR1 <- GO.MF.DR3vsDR1[GO.MF.DR3vsDR1$FDRunder<=0.05,]
write.table(enriched.under.GO.MF.DR3vsDR1, file = 'GO.molecular.function_under_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)

##### ANALISIS DE RUTAS KEGG #########

# creo FDRunder y FDRover
kegg_DE_DR3vsDR1$FDRunder <- p.adjust(kegg_DE_DR3vsDR1$under_represented_pvalue, n=nrow(kegg_DE_DR3vsDR1))
kegg_DE_DR3vsDR1$FDRover <- p.adjust(kegg_DE_DR3vsDR1$over_represented_pvalue, n=nrow(kegg_DE_DR3vsDR1))

kegg_path_DR3vsDR1<-mapPathwayToName("hsa")
idx_DR3vsDR1<-match(kegg_DE_DR3vsDR1$category,kegg_path_DR3vsDR1$path)
pathway_name_DR3vsDR1<-kegg_path_DR3vsDR1[idx_DR3vsDR1,2]
data_DR3vsDR1<-cbind(kegg_DE_DR3vsDR1,pathway_name_DR3vsDR1)

# Visualizo el resultado
head(data_DR3vsDR1)

write.table(data_DR3vsDR1, file = 'KEGG.PATHWAYS_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)

# guardo los que tienen un FDRover <= 0.05
enriched.KEGG.DR3vsDR1 <- data_DR3vsDR1[data_DR3vsDR1$FDRover<=0.05,]

# Genes under-expressed
under.KEGG.DR3vsDR1 <- data_DR3vsDR1[data_DR3vsDR1$FDRunder<=0.05,]
write.table(under.KEGG.DR3vsDR1, file = 'KEGG.PATHWAYS_under_DR3vsDR1.tsv', sep = "\t", row.names = FALSE)