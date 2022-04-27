library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb =TxDb.Mmusculus.UCSC.mm10.knownGene
library(org.Mm.eg.db)
library(clusterProfiler)


##  BMDM 

getwd()
setwd("~/Desktop/Bioinfo_analysis/Dawei_dualCHIP/BMDM_Spt16_GSE117333//")
BMDM_Spt16=readPeakFile("./idr_reproducibility/idr.conservative_peak.regionPeak.gz")
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_BMDM_Spt16 <- getTagMatrix(BMDM_Spt16, windows=promoter)


BMDM_Spt16_peakAnno <- annotatePeak(BMDM_Spt16, tssRegion=c(-3000, 3000),
                                    TxDb=txdb, annoDb="org.Mm.eg.db")
BMDM_Spt16_peakAnno_out=as.data.frame(BMDM_Spt16_peakAnno)

plotAnnoPie(BMDM_Spt16_peakAnno)
upsetplot(BMDM_Spt16_peakAnno)


setwd("~/Desktop/Bioinfo_analysis/Dawei_dualCHIP/BMDM_BRD4_GSE151015//")

BMDM_BRD4=readPeakFile("./idr_reproducibility/idr.conservative_peak.regionPeak.gz")
BMDM_BRD4_peakAnno <- annotatePeak(BMDM_BRD4, tssRegion=c(-3000, 3000),
                                   TxDb=txdb,annoDb="org.Mm.eg.db")
BMDM_BRD4_peakAnno_out=as.data.frame(BMDM_BRD4_peakAnno)

plotAnnoPie(BMDM_BRD4_peakAnno)
upsetplot(BMDM_BRD4_peakAnno)

library("ggvenn")

A <-list('BMDM_BRD4'=BMDM_BRD4_peakAnno_out$SYMBOL,'BMDM_Spt16'=BMDM_Spt16_peakAnno_out$SYMBOL)

# create venn diagram and display all the sets
ggvenn(A,set_name_size = 4,text_size =3,)

Mouse_IRG=read.csv("../IRG_mouse.csv")

BMDM_BRD4_IRG=intersect(BMDM_BRD4_peakAnno_out$SYMBOL,Mouse_IRG$x)
BMDM_Spt16_IRG=intersect(BMDM_Spt16_peakAnno_out$SYMBOL,Mouse_IRG$x)

D=list("BMDM_BRD4_IRG"=BMDM_BRD4_IRG,"BMDM_Spt16_IRG"=BMDM_Spt16_IRG)
ggvenn(D,set_name_size = 3,text_size =3)

getwd()
setwd("../ENCODE_results/")
BMDM_Spt16_Promoter=subset(BMDM_Spt16_peakAnno_out,BMDM_Spt16_peakAnno_out$annotation=="Promoter (<=1kb)"|BMDM_Spt16_peakAnno_out$annotation=="Promoter (1-2kb)"|BMDM_Spt16_peakAnno_out$annotation=="Promoter (2-3kb)")
write.csv(BMDM_Spt16_Promoter,"BMDM_Spt16_Promoter.csv")

BMDM_Brd4_Promoter=subset(BMDM_BRD4_peakAnno_out,BMDM_BRD4_peakAnno_out$annotation=="Promoter (<=1kb)"|BMDM_BRD4_peakAnno_out$annotation=="Promoter (1-2kb)"|BMDM_BRD4_peakAnno_out$annotation=="Promoter (2-3kb)")
write.csv(BMDM_Brd4_Promoter,"BMDM_Brd4_Promoter.csv")


BMDM_BRD4_Promoter_IRG=intersect(BMDM_Brd4_Promoter$SYMBOL,Mouse_IRG$x)
BMDM_Spt16_Promoter_IRG=intersect(BMDM_Spt16_Promoter$SYMBOL,Mouse_IRG$x)

D=list("BMDM_BRD4_Promoter_IRG"=BMDM_BRD4_Promoter_IRG,"BMDM_Spt16_Promoter_IRG"=BMDM_Spt16_Promoter_IRG)
ggvenn(D,set_name_size = 2,text_size =3)

E=list("BMDM_Brd4_Promoter"=BMDM_Brd4_Promoter$SYMBOL,"BMDM_Spt16_Promoter"=BMDM_Spt16_Promoter$SYMBOL)
ggvenn(E,set_name_size = 2,text_size =3)


#################### C2C12

setwd("~/Desktop/Bioinfo_analysis/Dawei_dualCHIP/C2C12_Spt16_GSE116169/")
C2C12_Spt16=readPeakFile("./idr_reproducibility/idr.conservative_peak.regionPeak.gz")

C2C12_Spt16_peakAnno <- annotatePeak(C2C12_Spt16, tssRegion=c(-3000, 3000),
                                     TxDb=txdb, annoDb="org.Mm.eg.db")
C2C12_Spt16_peakAnno_out=as.data.frame(C2C12_Spt16_peakAnno)

plotAnnoPie(C2C12_Spt16_peakAnno)
upsetplot(C2C12_Spt16_peakAnno)


setwd("~/Desktop/Bioinfo_analysis/Dawei_dualCHIP/C2C12_BRD4_GSE99101/")

C2C12_BRD4=readPeakFile("./idr_reproducibility/idr.conservative_peak.regionPeak.gz")
C2C12_BRD4_peakAnno <- annotatePeak(C2C12_BRD4, tssRegion=c(-3000, 3000),
                                    TxDb=txdb,annoDb="org.Mm.eg.db")
C2C12_BRD4_peakAnno_out=as.data.frame(C2C12_BRD4_peakAnno)


plotAnnoPie(C2C12_BRD4_peakAnno)

upsetplot(C2C12_BRD4_peakAnno)



C2C12_BRD4_Cistrome=readPeakFile("./Cistrome_83495_peaks.bed")

C2C12_BRD4_peakAnno_Cistrome <- annotatePeak(C2C12_BRD4_Cistrome, tssRegion=c(-3000, 3000),
                                             TxDb=txdb,annoDb="org.Mm.eg.db")

C2C12_BRD4_peakAnno_Cistrome_out=as.data.frame(C2C12_BRD4_peakAnno_Cistrome)

B=list('C2C12_BRD4'=C2C12_BRD4_peakAnno_out$SYMBOL,'C2C12_Spt16'=C2C12_Spt16_peakAnno_out$SYMBOL)

ggvenn(B,set_name_size = 4,text_size =3,)


C=list("C2C12"=C2C12_Spt16_peakAnno_out$SYMBOL,"BMDM"=BMDM_Spt16_peakAnno_out$SYMBOL)

ggvenn(C,set_name_size = 4,text_size =3,)



C2C12_BRD4_IRG=intersect(C2C12_BRD4_peakAnno_out$SYMBOL,Mouse_IRG$x)
C2C12_Spt16_IRG=intersect(C2C12_Spt16_peakAnno_out$SYMBOL,Mouse_IRG$x)

E=list("C2C12_BRD4_IRG"=C2C12_BRD4_IRG,"C2C12_Spt16_IRG"=C2C12_Spt16_IRG)
ggvenn(E,set_name_size = 3,text_size =3)


C2C12_Spt16_Promoter=subset(C2C12_Spt16_peakAnno_out,C2C12_Spt16_peakAnno_out$annotation=="Promoter (<=1kb)"|
                              C2C12_Spt16_peakAnno_out$annotation=="Promoter (1-2kb)"|
                              C2C12_Spt16_peakAnno_out$annotation=="Promoter (2-3kb)")


C2C12_BRD4_Promoter=subset(C2C12_BRD4_peakAnno_out,C2C12_BRD4_peakAnno_out$annotation=="Promoter (<=1kb)"|C2C12_BRD4_peakAnno_out$annotation=="Promoter (1-2kb)"|
                             C2C12_BRD4_peakAnno_out$annotation=="Promoter (2-3kb)")


C2C12_BRD4_Promoter_IRG=intersect(C2C12_BRD4_Promoter$SYMBOL,Mouse_IRG$x)
C2C12_Spt16_Promoter_IRG=intersect(C2C12_Spt16_Promoter$SYMBOL,Mouse_IRG$x)

D=list("C2C12_BRD4_Promoter_IRG"=C2C12_BRD4_Promoter_IRG,"C2C12_Spt16_Promoter_IRG"=C2C12_Spt16_Promoter_IRG)
ggvenn(D,set_name_size = 2,text_size =3)

E=list("C2C12_BRD4_Promoter"=C2C12_BRD4_Promoter$SYMBOL,"C2C12_Spt16_Promoter"=C2C12_Spt16_Promoter$SYMBOL)
ggvenn(E,set_name_size = 2,text_size =3)




############3 2022.04.06 reanalysis

gene=unique(BMDM_BRD4_peakAnno_out$SYMBOL)

