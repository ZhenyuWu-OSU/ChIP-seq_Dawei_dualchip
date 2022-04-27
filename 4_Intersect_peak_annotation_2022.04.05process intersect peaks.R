library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb =TxDb.Mmusculus.UCSC.mm10.knownGene
library(org.Mm.eg.db)
library(clusterProfiler)


getwd()
setwd("~/Desktop/Bioinfo_analysis/Dawei_dualCHIP/BEDtool intersect//")

intersections=readPeakFile("./BEDTOOL_intersect.bed.txt")

intersections_peakAnno <- annotatePeak(intersections, tssRegion=c(-3000, 3000),
                                   TxDb=txdb,annoDb="org.Mm.eg.db")
intersections_peakAnno_out=as.data.frame(intersections_peakAnno)

plotAnnoPie(intersections_peakAnno)
upsetplot(intersections_peakAnno)


Mouse_IRG=read.csv("../IRG_mouse.csv")

INTERSECT_IRG=subset(intersections_peakAnno_out,intersections_peakAnno_out$SYMBOL %in% Mouse_IRG$x)

unique_IRG=unique(INTERSECT_IRG$SYMBOL)

library("ggvenn")

intersections_2=readPeakFile("./BEDTOOL_intersect_reverse.txt")

intersections2_peakAnno <- annotatePeak(intersections_2, tssRegion=c(-3000, 3000),
                                       TxDb=txdb,annoDb="org.Mm.eg.db")
intersections_peakAnno2_out=as.data.frame(intersections2_peakAnno)



GALA=readPeakFile("./Galaxy4_interesect.bed")
GALA_ANNO= annotatePeak(GALA, tssRegion=c(-3000, 3000),
                        TxDb=txdb,annoDb="org.Mm.eg.db")

GALA_ANNO_out=as.data.frame(GALA_ANNO)

write.csv(GALA_ANNO_out,"peak_interesct.csv")

?plotAnnoPie
plotAnnoPie(GALA_ANNO,ndigit = 1,cex = 1.1)

# install.packages("UpSetR")
library("UpSetR")
UpSetR::upset(GALA_ANNO_out)
?upsetplot()
help(upsetplot)
upsetplot(GALA_ANNO)

unqie=unique(GALA_ANNO_out$SYMBOL)

GALA_ANNO_out_IRG=subset(GALA_ANNO_out,GALA_ANNO_out$SYMBOL %in% Mouse_IRG$x)
write.csv(GALA_ANNO_out_IRG,"peak_interesct_IRG.csv")

Intersect_Promoter=subset(GALA_ANNO_out,GALA_ANNO_out$annotation=="Promoter (<=1kb)"|GALA_ANNO_out$annotation=="Promoter (1-2kb)"|GALA_ANNO_out$annotation=="Promoter (2-3kb)")
write.csv(Intersect_Promoter,"peak_interesct_Promoter.csv")
unqie=unique(Intersect_Promoter_IRG$SYMBOL)


Intersect_Promoter_IRG=subset(Intersect_Promoter,Intersect_Promoter$SYMBOL %in% Mouse_IRG$x)
write.csv(Intersect_Promoter_IRG,"peak_interesct_Promoter_IRG.csv")






