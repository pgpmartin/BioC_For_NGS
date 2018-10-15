## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(kableExtra)

## ----warning=F,message=F-------------------------------------------------
library(Biostrings)

## ----cache=T-------------------------------------------------------------
dm3_upstream_filepath <- system.file("extdata","dm3_upstream2000.fa.gz",
                                    package="Biostrings")
dm3_upstream <- readDNAStringSet(dm3_upstream_filepath)
dm3_upstream
dm3_upstream[[5]]

## ---- cache=TRUE, warning=F, message=F-----------------------------------
library(Rsamtools)
fl <- system.file("extdata", "ce2dict1.fa", package="Rsamtools",
                  mustWork=TRUE)
fa <- open(FaFile(fl))
seqinfo(fa)

getSeq(fa, GRanges("pattern05:11-20"))

## ----warning=F,message=F-------------------------------------------------
library(BSgenome.Dmelanogaster.UCSC.dm3)

## ----cache=T-------------------------------------------------------------
names(Dmelanogaster)[1:5]
Dmelanogaster$chr2L

## ----eval=F--------------------------------------------------------------
## library(BSgenome.Dmelanogaster.UCSC.dm3.masked)

## ----eval=F--------------------------------------------------------------
## library(BSgenome)
## ?available.SNPs

## ----cache=T-------------------------------------------------------------
matchPattern("KATCGATA", dm3_upstream[[592]], fixed=FALSE)

## ----eval=F--------------------------------------------------------------
## vmatchPattern('TATCGATA', Dmelanogaster)

## ----echo=F,warning=F,message=F------------------------------------------
library(MotifDb);library(seqLogo)

## ----cache=T-------------------------------------------------------------
EcRMotif <- MotifDb::query(MotifDb,"EcR")[1]

## ----echo=F,cache=T,fig.height=4,fig.width=6,fig.align='center'----------
seqLogo::seqLogo(reverseComplement(EcRMotif[[1]]))

## ------------------------------------------------------------------------
EcRpfm <- apply(reverseComplement(EcRMotif[[1]]) *
                    as.integer(mcols(EcRMotif)$sequenceCount),
                2, as.integer)
rownames(EcRpfm) <- rownames(EcRMotif[[1]])
EcRpwm <- PWM(EcRpfm)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(EcRpwm,
             digits = rep(2,8),
             caption="PWM for EcR motif") %>%
    kable_styling(bootstrap_options = c("striped","condensed"),
                  full_width = FALSE,
                  font_size=14, position="center")

## ----warning=F,cache=T---------------------------------------------------
EcRHits <- matchPWM(EcRpwm, Dmelanogaster$chr4)
length(EcRHits)
EcRHits[1:2]

## ---- cache=T------------------------------------------------------------
pairwiseAlignment('CTTGCAGTGGTGTATTCATAC',
                  dm3_upstream[[1]],
                  type='global-local')

## ------------------------------------------------------------------------
stringDist(c("lazy", "HaZy", "crAzY"))
stringDist(c("lazy", "HaZy", "crAzY"), ignoreCase = TRUE)

## ------------------------------------------------------------------------
indx <- DNAStringSet( c("ATCACG", "CGATGT", "TTAGGC", "TGACCA",
                        "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA",
                        "GATCAG", "TAGCTT", "GGCTAC", "CTTGTA",
                        "CGGCTA", "TCCGCG", "ATGTCA", "AGCGAT"))

## ----echo=FALSE, warning=F,message=F-------------------------------------
library(ShortRead)

## ----cache=T-------------------------------------------------------------
fq1_path <- system.file(package="ShortRead","extdata","E-MTAB-1147",
                        "ERR127302_1_subset.fastq.gz")
myFastq <- readFastq(fq1_path)

## ------------------------------------------------------------------------
myFastq
myFastq[1:3]

## ------------------------------------------------------------------------
head(sread(myFastq), 2)
head(quality(myFastq), 2)
head(id(myFastq), 2)
encoding(quality(myFastq))[seq(1,51,by=2)]
alphabet(sread(myFastq))[1:4]

## ---- cache=T------------------------------------------------------------
nr_myFastq <- 0
strm <- FastqStreamer(fq1_path,1000)
repeat {
 ## Get FASTQ chunk:
  fq <- yield(strm)
  if (length(fq) == 0)
   break
 ## Do something on the chunk:
  nr_myFastq <- nr_myFastq + length(fq)
}
close(strm) #close the connection
nr_myFastq

## ----cache=T-------------------------------------------------------------
rqcResultSet <- rqcQA(fq1_path, sample=TRUE)

## ----cache=T, warning=F, out.width="500px", out.height="350px", fig.align='center'----
rqcCycleQualityPlot(rqcResultSet)

## ----cache=T, warning=F, out.width="500px", out.height="350px", fig.align='center'----
rqcCycleBaseCallsLinePlot(rqcResultSet)

## ----cache=T-------------------------------------------------------------
max1N <- nFilter(threshold=1L) #No 'Ns' in the reads
goodq <- srFilter(function(x){apply(as(quality(x),"matrix"),
                            1,median,na.rm=T)>=30},
                 name="MedianQualityAbove30")
myFilter <- compose(max1N,goodq) #combine filters

## ----cache=T-------------------------------------------------------------
FilterAndTrim <- function(fl,destination=sprintf("%s_filtered",fl))
{
  stream <- FastqStreamer(fl) ## open input stream
  on.exit(close(stream))
  repeat {
    ###get fastq chunk  
    fq <- yield(stream)
    if (length(fq)==0)
      break
    ###TRIM first 4 and last 2 bases
    fq <- narrow(fq,start=5,end=70)
    ###FILTER
    fq <- fq[myFilter(fq)]
    ###write filtered fastq
    writeFastq(fq, destination, mode="a")
  }
}

## ----eval=F--------------------------------------------------------------
## FilterAndTrim(fqFiles[1],
##               destination=file.path(getwd(),"FilteredFastq.fq"))

## ------------------------------------------------------------------------
library(Rsamtools);library(GenomicAlignments)
library(pasillaBamSubset)
sr <- untreated1_chr4() #single-end
pr <- untreated3_chr4() #paired-end

## ----eval=F--------------------------------------------------------------
## indexBam(sr_bamFile)

## ----cache=T-------------------------------------------------------------
which=GRanges(seqnames="chr4",
              ranges=IRanges(c(75000,1190000),
                             c(85000,1203000)),
              strand="*")
what = c("rname","strand","pos","qwidth","seq")
flag=scanBamFlag(isDuplicate=FALSE)
param=ScanBamParam(which=which,what=what,flag=flag)

## ----cache=T-------------------------------------------------------------
srbam <- readGAlignments(sr,param=param)
srbam[1:2]

## ----warning=F, cache=T--------------------------------------------------
prbam <- readGAlignmentPairs(pr)
prbam[1:2]

## ----cache=T-------------------------------------------------------------
eg = IRanges(start = c(1, 10, 20),
              end = c(4, 10, 19),
              names = c("A", "B", "C"))
eg

## ----cache=T-------------------------------------------------------------
set.seed(123) #For reproducibility
start = floor(runif(10000, 1, 1000))
width = floor(runif(10000, 0, 100))
ir = IRanges(start, width=width)
ir

## ----cache=T-------------------------------------------------------------
length(ir)
width(ir[1:4])
names(eg)

## ------------------------------------------------------------------------
mid(ir[1:4])
successiveIRanges(width=rep(10,3),gap=10)
tile(ir[1:2],n=2)

## ------------------------------------------------------------------------
irl  <- split(ir,width(ir)) # an IRangesList
irl[[1]][1:3]
length(irl)
head(elementNROWS(irl))

## ------------------------------------------------------------------------
start(irl)[1:2]
log(start(irl)[1:2])

## ------------------------------------------------------------------------
library(GenomicRanges)

## ----echo=FALSE----------------------------------------------------------
options(showHeadLines=3)
options(showTailLines=3)

## ----warning=F, message=F------------------------------------------------
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
(gr <- exons(txdb))

## ----cache=T-------------------------------------------------------------
(grl <- exonsBy(txdb,by="gene"))

## ----echo=F, eval=F------------------------------------------------------
## #plotRange function
## plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
## col = "black", sep = 0.5, cextitle=1, cexaxis=1, xaxis=T,...)
## {
## height <- 1
## if (is(xlim, "Ranges"))
## xlim <- c(min(start(xlim)), max(end(xlim)))
## bins <- disjointBins(IRanges(start(x), end(x) + 1))
## plot.new()
## plot.window(xlim, c(0, max(bins)*(height + sep)))
## ybottom <- bins * (sep + height) - height
## rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
## if (xaxis)
## (axis(1,cex.axis=cexaxis,padj=1))
## title(main,cex.main=cextitle)
## }

## ----echo=F,eval=F-------------------------------------------------------
## png("./figs/RangesOperations_intraRange.png",
##      width = 768, height = 768,res=72)
## par(mfrow=c(5,1))
## plotRanges((irs=IRanges(start=25,end=55)),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5,col=rep(c("blue","red"),each=3))
## plotRanges(shift(irs,40),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## plotRanges(resize(irs,20,fix="start"),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## plotRanges(promoters(irs,upstream=20,downstream=1),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## plotRanges(flank(irs,20,start=F),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## dev.off()
## 

## ---- echo=F, eval=F-----------------------------------------------------
## #Create irs
## irs = IRanges(start=c(1,9,25,46,53,87),
##               end=c(50,33,44,60,65,96))
## #Plot irs
## png("./figs/RangesOperations.png",
##      width = 850, height = 768,res=72)
## par(mfrow=c(5,1))
## plotRanges(irs,xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## plotRanges(reduce(irs),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## plotRanges(gaps(irs),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## plotRanges(disjoin(irs),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
## plot(1:100,c(coverage(irs),rep(0,4)),type="l",axes=F,xlab="",ylab="",lwd=3)
## title(main="coverage(irs)",cex.main=2)
## axis(side=2,lwd=2,cex.axis=2,at=0:3,labels=0:3)
## axis(1,lwd=2,cex.axis=2,padj=1)
## dev.off()

## ----echo=F,eval=F-------------------------------------------------------
## png("./figs/RangesOperations_setops.png",
##      width = 768, height = 768,res=72)
## par(mfrow=c(5,1))
## plotRanges(irs,xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5,col=rep(c("blue","red"),each=3))
## plotRanges(union(irs[1:3],irs[4:6]),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5,main="union")
## plotRanges(intersect(irs[1:3],irs[4:6]),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5,main="intersect")
## plotRanges(setdiff(irs[1:3],irs[4:6]),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5,col="blue",main="setdiff(blue,red)")
## plotRanges(setdiff(irs[4:6],irs[1:3]),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5,col="red",main="setdiff(red,blue)")
## dev.off()
## 

## ----cache=T-------------------------------------------------------------
Dmg <- genes(txdb) 
Dmt <- transcriptsBy(txdb,by="gene")

## ----cache=T-------------------------------------------------------------
Dm_tss <- unlist(reduce(promoters(Dmt,up=0,down=1),min.gap=0L))

## ----cache=T,warning=F---------------------------------------------------
mean(countOverlaps(Dm_tss,Dmg+500) > 1) #!strand-aware
mean(countOverlaps(Dm_tss,Dmg+500,ignore.strand=T) > 1)

## ----cache=T, warning=F--------------------------------------------------
fov <- findOverlaps(Dm_tss,Dmg+500,ignore.strand=T) ; fov[1:3]

## ----cache=T-------------------------------------------------------------
Dmg[subjectHits(fov)[queryHits(fov)==1]]

## ----cache=T-------------------------------------------------------------
length(subsetByOverlaps(srbam,Dmg["FBgn0002521"]))

## ----cache=T-------------------------------------------------------------
length(srbam[overlapsAny(srbam,grl[["FBgn0002521"]])])

## ----cache=T-------------------------------------------------------------
ctex <- summarizeOverlaps(features = grl[seqnames(Dmg)=="chr4"],
                             reads = srbam,
                              mode = Union)

## ----cache=T-------------------------------------------------------------
head(assays(ctex)$counts)

## ----cache=T-------------------------------------------------------------
srbam <- readGAlignments(sr)
(covr <- coverage(srbam))

## ----cache=T-------------------------------------------------------------
gn4 <- Dmg[seqnames(Dmg)=="chr4"]

## ----cache=T-------------------------------------------------------------
profgn4 <- covr[gn4]
profgn4[strand(gn4)=="-"] <- lapply(profgn4[strand(gn4)=="-"],rev)
names(profgn4) <- names(gn4) ; profgn4[1:2]

## ----cache=T-------------------------------------------------------------
profgn4 <- profgn4[elementNROWS(profgn4)>=1000]
profgn4 <- as(lapply(profgn4,window,1,1000),"RleList")
mat1kb <- matrix(as.numeric(unlist(profgn4, use.names=F)),
                 nrow=length(profgn4), byrow=T,
                 dimnames=list(names(profgn4),NULL))
mat1kb <- mat1kb[rowSums(mat1kb)>0,]

## ----cache=T-------------------------------------------------------------
df1Kb <- data.frame(Coordinate=1:1000,
                    Coverage=apply(mat1kb,2,mean,na.rm=T,trim=0.03))
ggplot(df1Kb,aes(x=Coordinate,y=Coverage))+
  geom_line()

## ----cache=T-------------------------------------------------------------
getSeq(Dmelanogaster,gn4[1:2])
Views(Dmelanogaster,gn4[1:2])

## ----echo=F,message=F,warning=F------------------------------------------
library(org.Dm.eg.db)

## ----message=F-----------------------------------------------------------
select(org.Dm.eg.db,
       keys=c('FBgn0015664','FBgn0015602'),keytype="FLYBASE",
       columns=c('SYMBOL','UNIGENE','ENTREZID','FLYBASECG'))

## ----eval=F--------------------------------------------------------------
## library(AnnotationHub)
## hub <- AnnotationHub()
## length(hub) # >43500 datasets
## unique(hub$dataprovider)
## head(unique(hub$species))
## head(unique(ah$rdataclass))

