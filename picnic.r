install.packages("devtools")
install.packages("BiocManager")
library("devtools")


BiocManager::install(c("TRONCO",
                       "TCGAbiolinks",
                       "maftools",
                       ))

install.packages("readr")

library(TCGAbiolinks)
library(maftools)
library(readr)
library(TRONCO)


source(functions.r)


BRCA.maf.file<-"input/data_mutations.txt"
BRCA.clinical.file<-"input/brca_tcga_clinical_data-2.tsv.txt"
BRCA.cna.file<-"input/data_cna.txt"
BRCA.mutex.file<-"input/ranked-groups.txt"

min_freq<-0.03
num_boot_iter<-100


#Apoptosis
CASP<-c("CASP8")
#Cell Cycle
CC<-c("CDKN1B", "RB1")
#CHROMATIN HISTONE MODIFIERS
CHM<-c("KMT2C","NCOR1")
#CHORMATION OTHER
CO<-c("CHD4","CTCF")
#CHROMATIN SWI/SNF COMPLEX
CSC<-c("ARID1A")
#GENOME INTEGRITY
TP53<-c("TP53","BRCA1")
#MAPK SIGNALING
MAPK<-c("KRAS","NF1",
        "CDH1","GPS2","MAP2K4","MAP3K1")
#PI3K SIGNALING
PI3K<-c("AKT1","PIK3CA","PIK3R1","PTEN",
        "PTPRD")
#PROTEIN HOMOESTASIUS/UBIQUITINATION
FBXW7<-c("FBXW7")
#RTK SIGNALING
ERBB2<-c("ERBB2")
#SPLICING
SF3B1<-c("SF3B1")
#TRANSCIPTION FACTOR
GATA<-c("CBFB","FOXA1","GATA3","RUNX1","TBX3")

pathway.genes <- c(
    CASP,
    CC,
    CHM,
    CO,
    CSC,
    TP53,
    MAPK,
    PI3K,
    FBXW7,
    ERBB2,
    SF3B1,
    GATA
    )
pathway.genes <- unique(pathway.genes)

pathway.names <- c(
    "CASP",
    "CC",
    "CHM",
    "CO",
    "CSC",
    "TP53",
    "MAPK",
    "PI3K",
    "FBXW7",
    "ERBB2",
    "SF3B1",
    "GATA"
    )

pathway.list <- list(
    CASP=CASP,
    CC=CC,
    CHM=CHM,
    CO=CO,
    CSC=CSC,
    TP53=TP53,
    MAPK=MAPK,
    PI3K=PI3K,
    FBXW7=FBXW7,
    ERBB2=ERBB2,
    SF3B1=SF3B1,
    GATA=GATA
    )

## colors for pathways in the plots
alteration.color <- 'dimgray'
pathways.color <- c(
    'darkslategray',
    'darkblue',
    'darkgreen',
    'darkmagenta',
    'firebrick4',
    'dodgerblue4',
    'darkorchid1',
    'darkorange4',
    'yellow',
    'darkred',
    'pink',
    'magenta'
    )

BRCA.maf<-maftools::read.maf(BRCA.maf.file)

#stampa un po' di informazioni
plotmafSummary(
  maf = BRCA.maf,rmOutlier = TRUE, addStat = 'median', dashboard = TRUE
)

#tutti i tipi di mutazioni in una sola
BRCA<- TRONCO::import.MAF(
    BRCA.maf.file,
    is.TCGA = TRUE,
    merge.mutation.types = TRUE,
    filter.fun = function(x){
        return(x['Hugo_Symbol'] %in% pathway.genes)
        }
    )

BRCA <- change.color(BRCA,
                     "Mutation",
                     "#D3AECC")
BRCA <- annotate.description(BRCA,
                             "BRCA somatic mutations")
oncoprint(BRCA)

clinical.data.brca<-TCGA.map.clinical.data(
  file = BRCA.clinical.file,
  column.samples = 'Patient.ID',
  column.map = 'Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code'
)
#levo le prime 3 righe
#pippo=pippo[-c(1,2,3,4),,drop=FALSE]
#rinomino la colonna
colnames(clinical.data.brca)[
    colnames(
        clinical.data.brca
        ) == "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code"
    ] ="pathologic_stage"

## bind together MAF and clinical data,
BRCA <- annotate.stages(BRCA,
                        clinical.data.brca,
                        match.TCGA.patients = TRUE)

## clearing the dataset
BRCA <- TCGA.remove.multiple.samples(BRCA)
BRCA <- TCGA.shorten.barcodes(BRCA)
BRCA <- annotate.stages(BRCA,
                        clinical.data.brca)

#print di mut+stage
oncoprint(BRCA)

#LOAD GISTIC FILE

#gisticfile
gistic <- read_delim(
    BRCA.cna.file,
    delim = "\t",
    escape_double = FALSE)
BRCA.gistic<-as.data.frame(gistic)

BRCA.gistic<-BRCA.gistic[BRCA.gistic$Hugo_Symbol %in% pathway.genes,]

BRCA.gistic<-t(BRCA.gistic)
colnames(BRCA.gistic)<-lapply(BRCA.gistic[1,],as.character)
BRCA.gistic<-BRCA.gistic[-1,]
BRCAGistic<-import.GISTIC(BRCA.gistic,trim = FALSE)
BRCAGistic<-annotate.description(BRCAGistic,label = "BRCA CNA data")
print(as.types(BRCAGistic))

oncoprint(BRCAGistic)

BRCAGistic<-delete.type(BRCAGistic,'Heterozygous Loss')
BRCAGistic<-delete.type(BRCAGistic,'Low-level Gain')
BRCAGistic<-rename.type(BRCAGistic, 'Homozygous Loss','Deletion')
BRCAGistic<-rename.type(BRCAGistic, 'High-level Gain','Amplification')
BRCAGistic$types['Deletion', ] <- '#8FBCBB'
BRCAGistic$types['Amplification', ] <- '#81a1c1'

oncoprint(BRCAGistic)


#associare GISTIC con clinical
BRCAGistic<-annotate.stages(
    BRCAGistic,
    clinical.data.brca,match.TCGA.patients = TRUE)

BRCAGistic<-TCGA.remove.multiple.samples(BRCAGistic)
BRCAGistic<-TCGA.shorten.barcodes(BRCAGistic)
BRCAGistic<-annotate.stages(BRCAGistic,clinical.data.brca)
oncoprint(BRCAGistic)

#intersect MAF and GISTIC
BRCA<-intersect.datasets(BRCAGistic,BRCA, intersect.genomes = FALSE)
BRCA<-trim(BRCA)
BRCA<-annotate.stages(BRCA,clinical.data.brca)
BRCA<-annotate.description(
  x=BRCA,
  label = "BRCA mutazioni somatiche e CNA"
  )
oncoprint(BRCA)

oncoprint(BRCA,
          legend.cex = .3,
          text.cex = 0.8,
          gene.annot = pathway.list,
          gene.annot.color = pathways.color)



#sottotipi
subtypes <- TCGAquery_subtype(tumor = "brca")

subtypes.Basal <- subtypes[subtypes$BRCA_Subtype_PAM50 %in% "Basal",]
subtypes.Her2 <- subtypes[subtypes$BRCA_Subtype_PAM50 %in% "Her2",]
subtypes.LumA <- subtypes[subtypes$BRCA_Subtype_PAM50 %in% "LumA",]
subtypes.LumB <- subtypes[subtypes$BRCA_Subtype_PAM50 %in% "LumB",]
subtypes.Normal <- subtypes[subtypes$BRCA_Subtype_PAM50 %in% "Normal",]

BRCA.Basal<-trim(samples.selection(BRCA, subtypes.Basal$patient))
BRCA.Her2<-trim(samples.selection(BRCA, subtypes.Her2$patient))
BRCA.LumA<-trim(samples.selection(BRCA, subtypes.LumA$patient))
BRCA.LumB<-trim(samples.selection(BRCA, subtypes.LumB$patient))
BRCA.Normal<-trim(samples.selection(BRCA, subtypes.Normal$patient))


BRCA.Basal<-annotate.description(BRCA.Basal, "BRCA Basal subtype")
BRCA.Her2<-annotate.description(BRCA.Her2, "BRCA Her2 subtype")
BRCA.LumA<-annotate.description(BRCA.LumA, "BRCA Luma subtype")
BRCA.LumB<-annotate.description(BRCA.LumB, "BRCA LumB subtype")
BRCA.Normal<-annotate.description(BRCA.Normal, "BRCA Normal subtype")

oncoprint(
    BRCA.LumA,
    gene.annot = pathway.list,gene.annot.color = pathways.color)
oncoprint(
    BRCA.LumB,
    gene.annot = pathway.list, gene.annot.color = pathways.color)
oncoprint(
    BRCA.Normal,
    gene.annot = pathway.list, gene.annot.color = pathways.color)
oncoprint(
    BRCA.Basal,
    gene.annot = pathway.list, gene.annot.color = pathways.color)
oncoprint(
    BRCA.Her2,
    gene.annot = pathway.list, gene.annot.color = pathways.color)

#mutual exclusivity
BRCA.mutex<-import.mutex.groups(
    BRCA.mutex.file
    )

#oncoprint dei primi due gruppi mutex
oncoprint(events.selection(BRCA,
                           filter.in.names = BRCA.mutex[[1]]),
          title = paste("BRCA - Mutex group 1"),
          legend.cex = .3,
          font.row = 4,
          cellheight = 10,
          ann.hits = FALSE,
          gene.annot = pathway.list,
          gene.annot.color = pathways.color)

oncoprint(events.selection(BRCA,
                           filter.in.names = BRCA.mutex[[2]]),
          title = paste("BRCA - Mutex group 2"),
          legend.cex = .3,
          font.row = 4,
          cellheight = 10,
          ann.hits = FALSE,
          gene.annot = pathway.list,
          gene.annot.color = pathways.color)

models<- list(BRCA,
              BRCA.Basal,
              BRCA.Her2,
              BRCA.LumA,
              BRCA.LumB,
              BRCA.Normal)

label.all<-'all subtypes'
label.basal<-'Basal'
label.her2<-'Her2'
label.luma<-'LumA'
label.lumb<-'LumB'
label.normal<-'Normal'

label.all.short<-'ALL'
label.basal.short<-'BAS'
label.her2.short<-'Her2'
label.luma.short<-'LumA'
label.lumb.short<-'LumB'
label.normal.short<-'NOR'

labels<-c(label.all,
          label.basal,
          label.her2,
          label.luma,
          label.lumb,
          label.normal)

labels.short<-c(label.all,
                label.basal.short,
                label.her2.short,
                label.luma.short,
                label.lumb.short,
                label.normal.short)


i <- 1
for (m in models) {
  
  # ricostruzione del modello
  # fatto in questo modo per permettere di analizzare le statistiche
  troncomodel <- model(m, labels[i], labels.short[i])
  
  i <- i + 1
}
