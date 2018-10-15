## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"------------------
library(knitr)
library(BiocStyle)
BiocStyle::latex()

## ----required_packages, echo=FALSE, message=FALSE,warning=FALSE----------
library(Biostrings)
library(BSgenome)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(BSgenome.Dmelanogaster.UCSC.dm3.masked)
library(knitr)
library(MotifDb,verbose=F)
library(seqLogo)
library(motifStack)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(ShortRead)
library(ggplot2)
library(Rqc)
library(pasillaBamSubset)
library(MMDiffBamSubset)
library(GenomeInfoDb)
library(Rsamtools)
library(GenomicAlignments)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(drosophila2.db)
library(drosophila2probe)
library(drosophila2cdf)
library(hom.Dm.inp.db)
library(GO.db)
library(rtracklayer)
library(AnnotationHub)
library(biomaRt)
library(SRAdb)
library(GEOquery)
library(Gviz)
#library(Rsubread)
#library(QuasR)
#library(ggbio)

## ----SaveSession,eval=FALSE----------------------------------------------
## save.image("MySession.RData")

## ----LoadSession,eval=FALSE----------------------------------------------
## load("MySession.RData")

## ----Sequence_analysis_title, eval=FALSE,echo=FALSE, message=FALSE-------
## ##---------------------------------------------------------##
## ## Manipulating strings and sequences in Bioconductor
## ##---------------------------------------------------------##

## ----SequenceContainers_title, eval=FALSE,echo=FALSE, message=FALSE------
## ##----------##
## ## Sequence containers and accessors
## ##--------- ##

## ----Biostring_lib,eval=FALSE,echo=FALSE---------------------------------
## require(Biostrings,warn.conflicts = F)

## ----dm3_upstream_path,cache=TRUE----------------------------------------
dm3_upstream_filepath = system.file("extdata",
                                    "dm3_upstream2000.fa.gz",
                                    package="Biostrings")

## ----dm3_upstream,cache=TRUE---------------------------------------------
dm3_upstream = readDNAStringSet(dm3_upstream_filepath)
dm3_upstream

## ----randomseq,cache=TRUE------------------------------------------------
randomSeq = DNAString(paste(sample(DNA_ALPHABET[1:4], 
                                   size=24,
                                   replace=TRUE), 
                            collapse=""))
randomSeq

## ----Dmelanogaster,eval=-1,cache=TRUE------------------------------------
library(BSgenome.Dmelanogaster.UCSC.dm3)
Dmelanogaster
names(Dmelanogaster)
Dmelanogaster$chr2L

## ----Access_dm3_upstream,cache=TRUE--------------------------------------
dm3_upstream[[5]]
toString(dm3_upstream[[5]][2:30])
subseq(dm3_upstream[[5]],start=2,end=30)
Views(dm3_upstream[[5]],start=c(1,11,21),end=c(10,20,30))

## ----FaFile_getSeq,eval=-1,cache=TRUE------------------------------------
library(Rsamtools)
indFaEx_path=system.file("extdata","ce2dict1.fa",package="Rsamtools")
indFaEx=FaFile(indFaEx_path)

getSeq(indFaEx,
       GRanges(c("pattern01:3-10",
                 "pattern04:10-24")))

## ----SequenceAnalysisMasks_title, eval=FALSE,echo=FALSE, message=FALSE----
## ##----------##
## ## Sequence analysis and masks
## ##--------- ##

## ----revcomp,cache=TRUE--------------------------------------------------
reverseComplement(dm3_upstream[[5]])

## ----alphabetFrequency,cache=TRUE----------------------------------------
alphabetFrequency(dm3_upstream[1:2],baseOnly=TRUE,as.prob=TRUE)

## ----letterFrequency,cache=TRUE------------------------------------------
letterFrequency(Dmelanogaster$chr2L,"CG",as.prob=TRUE)

## ----Dmel_genome_masked,eval=FALSE,echo=TRUE-----------------------------
## library("BSgenome.Dmelanogaster.UCSC.dm3.masked")

## ----Mask_activate,cache=TRUE--------------------------------------------
maskedgenome <- BSgenome.Dmelanogaster.UCSC.dm3.masked #A MaskedBSgenome object
chrU <- maskedgenome$chrU #A MaskedDNAString object
active(masks(chrU)) #Only some masks are active
active(masks(chrU)) <- TRUE #turn on all masks
chrUmask=injectHardMask(chrU) #Replaces the masked nucleotides by "+"
as(chrU,"XStringViews") #Get the unmasked regions

## ----Mask_views,cache=TRUE-----------------------------------------------
toString(Views(Dmelanogaster$chrU,start=1714848, width=12))
toString(Views(chrU,start=1714848, width=12))
toString(Views(chrUmask,start=1714848, width=12))

## ----Mask_matchPattern,cache=TRUE----------------------------------------
length(matchPattern('GAGAGAGAGAGA',maskedgenome$chrU))
length(matchPattern('GAGAGAGAGAGA',chrU))
length(matchPattern('GAGAGAGAGAGA',chrUmask))

active(masks(chrU))=F #deactivate all masks
length(matchPattern('GAGAGAGAGAGA',chrU))

active(masks(chrU))['RM']=T #activate only RepeatMasker
length(matchPattern('GAGAGAGAGAGA',chrU))

## ----PatternMatching_Alignment_title, eval=FALSE,echo=FALSE, message=FALSE----
## ##----------##
## ## Pattern matching and sequence alignment
## ##--------- ##

## ----MotifDb_lib,eval=TRUE,echo=FALSE------------------------------------
require(MotifDb)

## ----EcRmotif_query,cache=TRUE-------------------------------------------
EcRMotifs=MotifDb::query(MotifDb,"EcR")
EcRMotifs
EcRMotifs[[1]]

## ----seqLogo_lib,eval=FALSE,echo=FALSE-----------------------------------
## require(seqLogo,warn.conflicts = F, quietly = T)

## ----SeqLogo_Ecr,fig.cap="Sequence logo for EcR motif",fig.align="center",fig.width=6, fig.height=3,fig.pos="!h",out.width="0.5\\textwidth"----
seqLogo::seqLogo(EcRMotifs[[1]])

## ----SeqLogo_Ecr_revcomp,fig.cap="Sequence logo for EcR motif reverse complement",fig.align="center",fig.width=6, fig.height=3,fig.pos="!h",out.width="0.5\\textwidth"----
seqLogo::seqLogo(reverseComplement(EcRMotifs[[1]]))

## ----SeqLogo_EcrUsp,fig.cap="Sequence logo for EcR:Usp heterodimer",fig.align="center",fig.width=6, fig.height=3,fig.pos="!h",out.width="0.6\\textwidth"----
seqLogo::seqLogo(EcRMotifs[[2]])

## ----motifStack_Ecr,fig.cap="Sequence logos for EcR motifs using motifStack",fig.align="center",fig.width=6, fig.height=6,fig.pos="!h",out.width="0.6\\textwidth"----
#Format the PFMs:
pfms <- as.list(EcRMotifs)
names(pfms) <- gsub("Dmelanogaster-", "", names(pfms))
pfms <- lapply(names(pfms),
               function(x) {new("pfm",
                                mat=pfms[[x]],
                                name = x)})

#Plot
motifStack(pfms, layout="tree")

## ----EcrJASP,cache=TRUE--------------------------------------------------
EcrJASP <- EcRMotifs[[2]]

## ----EcRJASP_2PWM, cache=TRUE--------------------------------------------
nseq <- as.integer(mcols(EcRMotifs[2])$sequenceCount)
ecrpfm <- apply(round(nseq * EcrJASP,0), 2, as.integer)
rownames(ecrpfm) <- rownames(EcrJASP)
EcrJASP <- PWM(ecrpfm)

## ----warnOFF,eval=TRUE,echo=FALSE----------------------------------------
options(warn=-1)

## ----matchPWM_EcR_on2L,cache=TRUE----------------------------------------
EcRJASP_2L=matchPWM(EcrJASP,
                    Dmelanogaster$chr2L,
                    min.score='90%')
EcRJASP_2L_rev=matchPWM(reverseComplement(EcrJASP),
                        Dmelanogaster$chr2L,
                        min.score='90%')
#One way to merge the results:
EcRJASP_2L_all <- Views(Dmelanogaster$chr2L,
                        union(ranges(EcRJASP_2L),
                              ranges(EcRJASP_2L_rev)))

## ----warnON,eval=TRUE,echo=FALSE-----------------------------------------
options(warn=0)

## ----matchPWM_EcR_all, cache=TRUE, warning=FALSE-------------------------
EcRJASP_all <- matchPWM(EcrJASP, Dmelanogaster, min.score="90%")
EcRJASP_all

## ----TFBSTools_JASPAR2018_lib, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE----
library(TFBSTools)
library(JASPAR2018)

## ----EcR2018_fromJASPAR, cache=TRUE--------------------------------------
EcR2018 <- TFBSTools::getMatrixByID(JASPAR2018, "MA0534")

## ----SeqLogo_Ecr2018,fig.cap="Sequence logo for EcR motif (JASPAR2018)",fig.align="center",fig.width=6, fig.height=3,fig.pos="!h",out.width="0.5\\textwidth"----
TFBSTools::seqLogo(toICM(EcR2018))

## ------------------------------------------------------------------------
#hESR1:
ESR1 <- TFBSTools::getMatrixByID(JASPAR2018, "MA0112")
PWMSimilarity(toPWM(EcR2018), toPWM(ESR1), method = "Pearson")
#dCTCF
CTCF <- TFBSTools::getMatrixByID(JASPAR2018, "MA0531")
PWMSimilarity(toPWM(EcR2018), toPWM(CTCF), method = "Pearson")
#random permutation:
PWMSimilarity(toPWM(EcR2018), toPWM(permuteMatrix(EcR2018)), method = "Pearson")

## ----Ecr2018_on2L, warning=FALSE-----------------------------------------
EcR2018_on2L <- searchSeq(toPWM(EcR2018), Dmelanogaster$chr2L, min.score="90%")

## ----EcR_IR1_cons,cache=TRUE---------------------------------------------
EcR_IR1_cons=consensusString(EcrJASP,ambiguityMap="N")
EcR_IR1_cons
EcR_IR1_cons=substring(EcR_IR1_cons,first=2)
EcR_IR1_cons

## ----EcR_Usp_cons,cache=TRUE---------------------------------------------
EcR_cons=substring(consensusString(reverseComplement(EcRMotifs[[1]]),
                                   ambiguityMap="N"),
                   first=2)
Usp_cons=substring(consensusString(MotifDb::query(MotifDb,"Usp")[[5]],
                                   ambiguityMap="N"),
                   first=1,
                   last=7)

## ----matchPattern_EcR_on2L,cache=TRUE------------------------------------
EcR_on_2L=matchPattern(EcR_cons,Dmelanogaster$chr2L)
EcR_on_2L_rev=matchPattern(reverseComplement(DNAString(EcR_cons)),
                           Dmelanogaster$chr2L)
EcR_on_2L_all=Views(Dmelanogaster$chr2L,union(ranges(EcR_on_2L),
                                              ranges(EcR_on_2L_rev)))

## ----Perfect_palindrome_on2L,cache=TRUE----------------------------------
EcR_on_2L_all[width(EcR_on_2L_all)!=7]

## ----vmatchPattern_EcR_dmUp,cache=TRUE-----------------------------------
EcR_on_up=vmatchPattern(EcR_cons, dm3_upstream)

## ----nmatch_per_seq,cache=TRUE-------------------------------------------
nmatch_per_seq = elementNROWS(EcR_on_up)
table(nmatch_per_seq)

## ----max_matches,cache=TRUE----------------------------------------------
i0=which.max(nmatch_per_seq)
Views(dm3_upstream[[i0]], EcR_on_up[[i0]])

## ----dm_matrices,cache=TRUE----------------------------------------------
dm_matrices = MotifDb::query(MotifDb,"dmelanogaster")

## ----dm_motifs,cache=TRUE------------------------------------------------
motif_ln = sapply(dm_matrices,ncol)
dm_matrices = dm_matrices[motif_ln==8]
dm_motifs=DNAStringSet(sapply(dm_matrices,consensusString,ambiguityMap="N"))

## ----matchPDict_on2L,cache=TRUE------------------------------------------
mot8_on_2L=matchPDict(dm_motifs,Dmelanogaster$chr2L,fixed=FALSE)
summary(elementNROWS(mot8_on_2L)) #Number of matches
head(unlist(mot8_on_2L)) #first 6 matches

## ----Frequent_motif,cache=TRUE-------------------------------------------
names(dm_motifs[which.max(elementNROWS(mot8_on_2L))])
toString(dm_motifs[[which.max(elementNROWS(mot8_on_2L))]])

## ----Rare_motif,cache=TRUE-----------------------------------------------
names(dm_motifs[which.min(elementNROWS(mot8_on_2L))])
toString(dm_motifs[[which.min(elementNROWS(mot8_on_2L))]])

## ----dm_mot8_dict,cache=TRUE---------------------------------------------
dm_mot8_dict=PDict(dm_motifs[sapply(dm_motifs,hasOnlyBaseLetters)])

## ----vcountPDict,cache=TRUE----------------------------------------------
mot8_count_upstream=vcountPDict(dm_mot8_dict,dm3_upstream)

## ----Number_of_motifs_found,cache=TRUE-----------------------------------
apply(mot8_count_upstream,1,sum)

## ----Number_of_motif1_per_seq,cache=TRUE---------------------------------
table(mot8_count_upstream[1,])

## ----nMot8_perUpSeq,cache=TRUE-------------------------------------------
nMot8_perSeq=apply(mot8_count_upstream,2,sum)
names(nMot8_perSeq)=names(dm3_upstream)

## ----fig_Mot8PerUpSeq,fig.align="center",fig.width=5, fig.height=5,fig.cap='Motifs per upstream sequence.',fig.pos="!h",out.width="0.6\\textwidth"----
nMot8_perSeq=nMot8_perSeq[nMot8_perSeq>=1]
plot(as.integer(names(table(nMot8_perSeq))),
  	as.integer(table(nMot8_perSeq)),
		pch=15,type="b",log="y",col="blue",
		xlab="Number of 8bp-long motifs",
		ylab="Number of upstream sequences")

## ----pairwiseAlignment_ex1,cache=TRUE------------------------------------
pairwiseAlignment(EcR_IR1_cons,
                  Dmelanogaster$chr2L[14067:14101],
                  type='global-local')

## ----pairwiseAlignment_ex2,cache=TRUE------------------------------------
paln_EcRUsp=pairwiseAlignment(c(EcR_cons,Usp_cons),
                              dm3_upstream[[1780]],
                              type='global-local')
paln_EcRUsp[1]
paln_EcRUsp[2]
Views(paln_EcRUsp)
Views(dm3_upstream[[1780]],start=772,end=785)

## ----pairwiseAlignment_ex3,cache=TRUE------------------------------------
pairwiseAlignment('GTGTCAATACGACAGCAATCTG',
                  'AGTGTGAATTACAGCAAATCTCTGTT',
                  type='global')
pairwiseAlignment('GTGTCAATACGACAGCAATCTG',
                  'AGTGTGAATTACAGCAAATCTCTGTT',
                  type='local')
pairwiseAlignment('GTGTCAATACGACAGCAATCTG',
                  'AGTGTGAATTACAGCAAATCTCTGTT',
                  type='global-local')

## ----pairwiseAlignment_ex4,cache=TRUE------------------------------------
pairwiseAlignment('GTGTCAATACGACAGCAATCTG',
                  'AGTGTGAATTACAGCAAATCTCTGTTCAATTTCTG',
                  type='global')
pairwiseAlignment('GTGTCAATACGACAGCAATCTG',
                  'AGTGTGAATTACAGCAAATCTCTGTTCAATTTCTG',
                  gapExtension=-6,type='global')
pairwiseAlignment('GTGTCAATACGACAGCAATCTG',
                  'AGTGTGAATTACAGCAAATCTCTGTTCAATTTCTG',
                  gapOpening=-60,type='global')

## ----pairwiseAlignment_utilsAndAccess,cache=TRUE-------------------------
paln=pairwiseAlignment(
               DNAStringSet(c('GTGTCAATACGACAGCAATCTG',
                              'TAAGGTCATAGTGT')),
               DNAString('TCGCCATAGGTCAATAGTGTGAATTACAGCAAATCTCTGTTCAATTTCTG'),
                       type='global-local')
pattern(paln)
subject(paln)
aligned(paln)
Biostrings::score(paln)
pid(paln) # percentage identity
compareStrings(paln) #symbolic representation of the alignment
nedit(paln) #Levenshtein edit distance
nmatch(paln) #also nmismatch(paln)
insertion(paln)[[1]] #also deletion(paln) 
Views(paln)
coverage(paln)

## ----stringDist_Function, cache=TRUE-------------------------------------
paln[1]
nedit(paln[1])
#Can also be obtained with:
stringDist(c('GTGTCAATACGACAGCAATCTG', 
             'GTGTGAATTACAGCAAATCTC'), 
           method = "levenshtein")

## ----pairwiseAlignment_protein,cache=TRUE--------------------------------
data(BLOSUM62)
paln <- pairwiseAlignment(AAString("ALAKHLYDSYIKSFPLTKAKARAILTGKTTDKS"),
                          AAString("AQQFNDIVCAMTQEDLEKFWKRCSRPFTAHM"),
                          substitutionMatrix = BLOSUM62)

## ----DECIPHER BrowseSeqs, eval=FALSE-------------------------------------
## library(DECIPHER)
## alnseqs <- c(aligned(pattern(paln)), aligned(subject(paln)))
## BrowseSeqs(alnseqs)

## ----GRanges_title, eval=FALSE,echo=FALSE, message=FALSE-----------------
## ##---------------------------------------------------------##
## ## Manipulating genomic ranges
## ##---------------------------------------------------------##

## ----IRangesIntro_title, eval=FALSE,echo=FALSE, message=FALSE------------
## ##----------##
## ## IRanges and accessors
## ##--------- ##

## ----IRanges_lib,eval=FALSE,echo=FALSE-----------------------------------
## require(IRanges,warn.conflicts = F, quietly = T)

## ----simple_IRanges,cache=TRUE-------------------------------------------
eg = IRanges(start = c(1, 10, 20),
              end = c(4, 10, 19),
              names = c("A", "B", "C"))
eg

## ----bigger_IRanges,cache=TRUE-------------------------------------------
set.seed(123) #For reproducibility
start = floor(runif(10000, 1, 1000))
end = start + floor(runif(10000, 0, 100))
ir = IRanges(start, end)
ir

## ----IRanges_accessors,cache=TRUE----------------------------------------
length(ir)
ir[1:4]
start(ir[1:4])
width(ir[1:4])
names(eg)

## ----IRanges_otherMethods,cache=TRUE-------------------------------------
c(ir[1:2],ir[5:6]) #combining
sort(ir[1:4])
rank(ir[1:4],ties="first")
mid(ir[1:4]) # midpoints
tile(ir[1:2],n=2) #returns an IRangesList (see below)
ir[[1]]
as.integer(ir[1]) #equivalent but works on multiple ranges
unlist(ir[1]) #also equivalent but names can be added
rep(ir[1:2],each=2)
isNormal(ir[1:4])
isNormal(sort(ir[1:4])) #see ?'Ranges-class' for Normality definition
isDisjoint(ir[1:4])
match(ir[1:4],ir[4:1]) #see ?'Ranges-comparison' for Ranges comparison methods
ir[1:4]>ir[4:1]

## ----otherIRanges_functions,cache=TRUE-----------------------------------
as(c(2:10,8,90:100),"IRanges") #from a vector of integers
successiveIRanges(width=rep(10,5),gap=10)
whichAsIRanges(c(19, 5, 0, 8, 5)>=5) #transforms a logical vector in IRanges

## ----IRangesList,cache=TRUE----------------------------------------------
irl=split(ir,width(ir)) # an IRangesList
irl[[1]]
start(irl)
head(elementNROWS(irl))

## ----IntraInterRange_ops_title, eval=FALSE,echo=FALSE, message=FALSE-----
## ##----------##
## ## Intra-and inter-range operations
## ##--------- ##

## ----intra-range_ops,cache=TRUE------------------------------------------
ir[1:2]
shift(ir[1:2],shift=10)
resize(ir[1:2],width=100,fix="start")
flank(ir[1:2],width=100,start=T)
narrow(ir[1:2],start=1,width=30) #here 'start' is relative
ir[1:2]+10
ir[1:2]-10
ir[1:2]+c(0,10)
ir[1:4]*-10 ; ir[1:4]*10 # acts like a centered zoom
ir[1:2]*c(1,2) #zoom second range by 2X

## ----plotRanges_function,cache=TRUE--------------------------------------
plotRanges <- function(x, 
                       xlim = x, 
                       main = deparse(substitute(x)),
                       col = "black", sep = 0.5, 
                       cextitle=1, cexaxis=1, xaxis=T,...) 
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  if (xaxis)
    (axis(1,cex.axis=cexaxis,padj=1))
  title(main,cex.main=cextitle)
}

## ----inter-ranges_ops,cache=TRUE-----------------------------------------
#select some ranges in [1:100]:
irs=ir[which(start(ir)<=100 & 
               end(ir)<=100)[c(3:4,8,14,17:18)]] 
irs
reduce(irs)
disjoin(irs)
gaps(irs)
coverage(irs)

## ----fig_inter-range_ops,fig.align="center",fig.width=4, fig.height=8,fig.cap='Inter-range operations.',fig.pos="!h",out.width="5in",out.height="7in"----
par(mfrow=c(5,1))
plotRanges(irs,xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
plotRanges(reduce(irs),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
plotRanges(disjoin(irs),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
plotRanges(gaps(irs),xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
plot(1:100,c(coverage(irs),rep(0,4)),type="l",axes=F,xlab="",ylab="",lwd=3)
title(main="coverage",cex.main=2)
axis(side=2,lwd=2,cex.axis=2,at=0:3,labels=0:3)
axis(1,lwd=2,cex.axis=2,padj=1)

## ----setops_IRanges_unionintersect,cache=TRUE----------------------------
union(irs[1:3],irs[4:6])
intersect(irs[1:3],irs[4:6])

## ----fig_IRanges_setops_unionintersect,fig.align="center",fig.width=4, fig.height=8, fig.cap='Union and intersect on IRanges.',fig.pos="!h",out.width="5in",out.height="5in"----
par(mfrow=c(4,1))
 plotRanges(irs[1:3], xlim=c(0,100), xaxis=T, 
            cextitle=2, cexaxis=1.5)
 plotRanges(irs[4:6], xlim=c(0,100), xaxis=T, 
            cextitle=2, cexaxis=1.5)
 plotRanges(union(irs[1:3], irs[4:6]), xlim=c(0,100), xaxis=T, 
            cextitle=2, cexaxis=1.5)
 plotRanges(intersect(irs[1:3], irs[4:6]), xlim=c(0,100), xaxis=T, 
            cextitle=2, cexaxis=1.5)

## ----setops_IRanges_setdiff,cache=TRUE-----------------------------------
setdiff(irs[1:3],irs[4:6])
setdiff(irs[4:6],irs[1:3])

## ----fig_IRanges_setops_setdiff,fig.align="center",fig.width=4, fig.height=8, fig.cap='Asymetric differences with setdiff on IRanges.',fig.pos="!h",out.width="5in",out.height="5in"----
par(mfrow=c(4,1))
 plotRanges(irs[1:3],xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
 plotRanges(irs[4:6],xlim=c(0,100),xaxis=T, cextitle=2, cexaxis=1.5)
 plotRanges(setdiff(irs[1:3],irs[4:6]),xlim=c(0,100),xaxis=T, 
            cextitle=2, cexaxis=1.5)
 plotRanges(setdiff(irs[4:6],irs[1:3]),xlim=c(0,100),xaxis=T, 
            cextitle=2, cexaxis=1.5)

## ----setops_IRanges_paralell,cache=TRUE----------------------------------
punion(irs[1:2],irs[4:5]) #element-wise (aka "parallel") union
pintersect(irs[1:2],irs[4:5])
psetdiff(irs[1:3],irs[4:6]) # asymmetric! difference
pgap(irs[1:3],irs[4:6])

## ----nearest_IRanges,cache=TRUE------------------------------------------
nearest(irs[4:6],irs[1:3])
distance(irs[4:6],irs[1:3])

## ----Rle_title, eval=FALSE,echo=FALSE, message=FALSE---------------------
## ##----------##
## ## Rle
## ##--------- ##

## ----Rle_objects,cache=TRUE----------------------------------------------
set.seed(123)
lambda = c(rep(0.001, 3500), #From IRanges vignette
           seq(0.001, 10, length = 500), 
            seq(10, 0.001, length = 500))
xRle=Rle(rpois(1e4,lambda))
yRle=Rle(rpois(1e4, lambda[c(251:length(lambda), 1:250)]))
xRle
yRle
as.vector(object.size(xRle)/object.size(as.vector(xRle))) #Gain of memory
head(runValue(xRle))
head(runLength(xRle))
head(start(xRle)) #starts of the runs
head(end(xRle)) #ends of the runs
nrun(xRle) #number of runs
findRun(as.integer(c(100,200,300,1200)),xRle)
coverage(irs)

## ----Rle_operations,cache=TRUE-------------------------------------------
xRle+yRle
xRle>0
xRle>yRle
max(xRle)
summary(xRle)
sqrt(xRle)
rev(xRle)
table(xRle)
union(xRle,yRle)
cor(xRle,yRle)

## ----Rle_runningWindow,cache=TRUE----------------------------------------
runmean(xRle,k=100) # See ?'Rle-runstat' for other examples
#same result, more flexible but much slower:
Rle(aggregate(xRle, start = 1:(length(xRle)-99), width = 100, FUN = mean))
runq(xRle,k=100,i=10) #10th smallest value in windows of 100

## ----Views_title, eval=FALSE,echo=FALSE, message=FALSE-------------------
## ##----------##
## ## Views
## ##--------- ##

## ----RleViews_create,cache=TRUE------------------------------------------
coverage(irs)
irs_views=Views(coverage(irs),start=c(-5,10,20,90),end=c(10,30,50,100))
irs_views #Views can be out of bound
try(irs_views[[1]]) #but can't be extracted
irs_views[[2]]
start(irs_views)

## ----RleViews_create2,cache=TRUE-----------------------------------------
Views(coverage(irs),irs[c(1,2,6)]) #use an IRanges to extract Views
Views(coverage(irs),coverage(irs)>=1) #or a logical Rle
slice(coverage(irs),3) #use slice
successiveViews(coverage(irs),width=rep(20,4)) #get successive Views

## ----XStringViews,cache=TRUE---------------------------------------------
dmup_views=Views(dm3_upstream[[1]],irs[1:2])
dmup_views
nchar(dmup_views)
toString(dmup_views) #see ?'XStringViews-class' for other methods

## ----RleViewsList,cache=TRUE---------------------------------------------
xyRleList=RleList(xRle,yRle)
xyRleList_views=Views(xyRleList,
                      IRangesList(IRanges(start=c(3700,4000),width=20),
                                  IRanges(yRle>17)+2))
xyRleList_views[[1]]
xyRleList_views[[2]]
width(xyRleList_views)

## ----Views_looping,cache=TRUE--------------------------------------------
viewMins(irs_views) #same as min(irs_views)
viewSums(irs_views) #same as sum(irs_views)
viewWhichMaxs(irs_views) #get the (first) coordinate of viewMaxs 
# (which.max also works)
viewRangeMins(irs_views) #get the (first) range of viewMins
viewApply(irs_views,sd)
viewMeans(xyRleList_views)

## ----GenomicRanges_title, eval=FALSE,echo=FALSE, message=FALSE-----------
## ##----------##
## ## GenomicRanges
## ##--------- ##

## ----GenomicRanges_lib,eval=FALSE,echo=FALSE-----------------------------
## require(GenomicRanges,warn.conflicts = F, quietly = T)

## ----GRanges_create,cache=TRUE-------------------------------------------
genes = GRanges(seqnames=c("chr2L", "chrX"),
                ranges=IRanges(start=c(7529, 18962306),
                               end =c(9484, 18962925),
                               names=c("FBgn0031208", "FBgn0085359")),
                strand=c("+", "-"),
                seqlengths=c(chr2L=23011544L, chrX=22422827L))

slotNames(genes)

mcols(genes) = DataFrame(EntrezId=c("33155", "2768869"),
                         Symbol=c("CG11023", "CG34330"))
genome(genes)="dm3" #see ?seqinfo for details
genes

## ----GRanges_accessors,cache=TRUE----------------------------------------
width(genes)
names(genes)
seqnames(genes)
strand(genes)
ranges(genes)
genes$Symbol
mcols(genes)
seqinfo(genes)
seqlevels(genes)

## ----GRangesList_objects,cache=TRUE--------------------------------------
gr1 = GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
              strand = "+", score = 5L, GC = 0.45)
gr2 = GRanges(seqnames = c("chr1", "chr1"),
              ranges = IRanges(c(7,13), width = 3),
              strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
grl = GRangesList("txA" = gr1, "txB" = gr2)
grl
length(grl)
elementNROWS(grl)
grl["txB","GC"]
unlist(grl)

## ----GenomicFeatures_lib,eval=FALSE,echo=FALSE---------------------------
## require(GenomicFeatures, quietly = T)

## ----TxDb_Dmel_lib,eval=FALSE,echo=FALSE---------------------------------
## require(TxDb.Dmelanogaster.UCSC.dm3.ensGene, quietly = T)

## ----TxDb_GRanges_extract,eval=TRUE,echo=TRUE,cache=TRUE-----------------
Dmg=genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene,single.strand.genes.only=T)
Dmg

## ----TxDb_GRangesList_extract,eval=TRUE,echo=TRUE,cache=TRUE-------------
Dmt=transcriptsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene,by="gene")
Dmt

## ----TxDb_extractByOverlap,cache=TRUE------------------------------------
exonsByOverlaps(TxDb.Dmelanogaster.UCSC.dm3.ensGene,genes)

## ----TxDb_convert2Seq,cache=TRUE-----------------------------------------
getSeq(BSgenome.Dmelanogaster.UCSC.dm3,Dmg[1:2])
Dmc=cdsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene,by="tx")
cds_seq=extractTranscriptSeqs(BSgenome.Dmelanogaster.UCSC.dm3,Dmc[1:2])
cds_seq
translate(cds_seq)

## ----GRanges_methods,cache=TRUE------------------------------------------
genes
genes2=Dmg[c(1:2,21:22,36:37)]
genes2
sort(genes2)
c(genes,genes2,ignore.mcols=T) #combine
intersect(genes,c(genes2,Dmg['FBgn0031208'])) #set operations
nearest(genes,genes2) #nearest-methods
nearest(genes,genes2,ignore.strand=T) #strand-aware by default
precede(genes,genes2,ignore.strand=T)
promoters(genes,upstream=200,downstream=1) #intra-range operations
reduce(Dmt[[2]]) #inter-range operations

## ----countOverlaps_example,eval=TRUE,echo=TRUE,message=FALSE,cache=TRUE,warning=FALSE----
Dm_tss=unlist(reduce(promoters(Dmt,up=0,down=1))) #get all TSS
cov_tss_g500=countOverlaps(Dm_tss,Dmg+500) #strand-aware!
table(cov_tss_g500)
sum(cov_tss_g500>1)
cov_tss_g500_bs=countOverlaps(Dm_tss,Dmg+500,ignore.strand=T) #both strands
sum(cov_tss_g500_bs>1)

## ----findOverlaps_example,eval=TRUE,echo=TRUE,message=FALSE,cache=TRUE,warning=FALSE----
fov_tss_g500_bs=findOverlaps(Dm_tss,Dmg+500,ignore.strand=T)
Dmg[c(1,1383)]

## ----randomreads2L,cache=TRUE--------------------------------------------
set.seed(0)
randomreads2L=GRanges(seqnames="chr2L",
                    ranges=IRanges(start=floor(runif(10000,5000,50000)),
                                   width=100),
                    strand="*")

## ----Dmg_2genes,cache=TRUE-----------------------------------------------
sort(Dmg)[1:2]

## ----subsetByOverlaps_example,cache=TRUE---------------------------------
subsetByOverlaps(randomreads2L,sort(Dmg)[1:2])

## ----overlapsAny_example,cache=TRUE--------------------------------------
randomreads2L[overlapsAny(randomreads2L,sort(Dmg)[1:2])]

## ----counting reads_example,cache=TRUE-----------------------------------
assays(summarizeOverlaps(sort(Dmg)[1:2],randomreads2L,mode="Union"))$counts

## ----RangedData_title, eval=FALSE,echo=FALSE, message=FALSE--------------
## ##----------##
## ## RangedData
## ##--------- ##

## ----Fastq_title, eval=FALSE,echo=FALSE, message=FALSE-------------------
## ##---------------------------------------------------------##
## ## Working with Fastq files
## ##---------------------------------------------------------##

## ----FASTQ_format_title, eval=FALSE,echo=FALSE, message=FALSE------------
## ##----------##
## ## FASTQ format
## ##--------- ##

## ----FASTQ_reading_title, eval=FALSE,echo=FALSE, message=FALSE-----------
## ##----------##
## ## Reading FASTQ files
## ##--------- ##

## ----ShortRead_lib,eval=FALSE,echo=FALSE---------------------------------
## require(ShortRead, quietly = T)

## ----ggplot2_lib,eval=FALSE,echo=FALSE-----------------------------------
## require(ggplot2, quietly = T)

## ----readFastq_example,eval=TRUE,echo=TRUE,cache=TRUE--------------------
fq1_path=system.file(package="ShortRead","extdata","E-MTAB-1147",
                    "ERR127302_1_subset.fastq.gz")
# A FASTQ file containing 20K reads
myFastq=readFastq(fq1_path)

## ----ShortReadQ_accessors,cache=TRUE-------------------------------------
myFastq
myFastq[1:5]
head(sread(myFastq),3)
head(quality(myFastq),3)
head(id(myFastq),3)
encoding(quality(myFastq))[2:50]
alphabet(sread(myFastq))

## ----FastqSampler,cache=TRUE---------------------------------------------
set.seed(123)
fqs1K = FastqSampler(fq1_path,1000)
reads_sample=yield(fqs1K)
close(fqs1K) #close connection
reads_sample

## ----FastqStreamer,cache=TRUE--------------------------------------------
nr_myFastq=0
strm <- FastqStreamer(fq1_path,1000)
repeat {
        fq <- yield(strm)
        if (length(fq) == 0)
          break
  ## Get FASTQ chunk
        nr_myFastq=nr_myFastq+length(fq)
  ## Do something on the chunk 
}
close(strm) #close the connection
nr_myFastq

## ----QAonFASTQ_title, eval=FALSE,echo=FALSE, message=FALSE---------------
## ##----------##
## ## Quality assessment on FASTQ files
## ##--------- ##

## ----ShortRead_qa2,cache=TRUE--------------------------------------------
fqPath = system.file(package="ShortRead", "extdata", "E-MTAB-1147")
fqFiles = dir(fqPath, pattern="fastq.gz", full=TRUE)
coll = QACollate(QAFastqSource(fqFiles), QAReadQuality(),
                  QAAdapterContamination(), QANucleotideUse(),
                  QAQualityUse(), QASequenceUse(),
                  QAFrequentSequence(n=10), QANucleotideByCycle(),
                  QAQualityByCycle())
qa2OnFastq = qa2(coll,BPPARAM=SerialParam(), verbose=FALSE)
## qa_report=report(qa2OnFastq) #generate the report
## browseURL(qa_report) #display in your browser
slotNames(qa2OnFastq)
names(qa2OnFastq)
slotNames(qa2OnFastq[["QANucleotideUse"]])
qa2OnFastq[["QANucleotideUse"]]

## ----ShortRead_qa,cache=TRUE,fig.align="center",fig.width=8, fig.height=5,fig.cap='Distribution of average base quality.',fig.pos="!h",out.width="0.8\\textwidth"----
qaOnFastq = qa(fqPath,"fastq.gz",BPPARAM=SerialParam()) #collect statistics
qaOnFastq[["readCounts"]]
readQuals <- qaOnFastq[["readQualityScore"]]
ggplot(readQuals,aes(x=quality, y=density, colour=lane))+
    geom_line(size=1)+
    scale_colour_discrete(name="Lane",
                          breaks=c("ERR127302_1_subset.fastq.gz",
                                   "ERR127302_2_subset.fastq.gz"),
                          labels=c("ERR127302_1","ERR127302_2")) +
    theme_bw()

## ----Rqc_lib,eval=FALSE,echo=FALSE---------------------------------------
## require(Rqc, quietly = T)

## ----Rqc_examples,fig.width=5, fig.height=5,fig.pos="h",out.width="0.49\\linewidth",out.height="6cm",fig.show='hold',fig.cap="Examples of QC plots using the Rqc package",cache=TRUE,warning=F----
rqcResultSet = rqcQA(fqFiles,sample=T)
rqcCycleQualityPlot(rqcResultSet[1])
rqcCycleBaseCallsLinePlot(rqcResultSet[2])

## ----Reads_fitering_title, eval=FALSE,echo=FALSE, message=FALSE----------
## ##----------##
## ## Reads filtering and trimming
## ##--------- ##

## ----ShortRead_create_filters,cache=TRUE---------------------------------
max1N=nFilter(threshold=1L)
#Remove reads with more than 1N
goodq = srFilter(function(x){
                            apply(as(quality(x),"matrix"),
                            1,median,na.rm=T)>=30
                            },
                 name="MedianQualityAbove30")
#Custom filter: Remove reads with median quality<30
myFilter=compose(max1N,goodq) #combine filters

## ----ShortRead_filterfunction,eval=TRUE,cache=TRUE-----------------------
FilterAndTrim = function(fl,destination=sprintf("%s_filtered",fl))
{
  stream = FastqStreamer(fl)
  on.exit(close(stream))
  ## open input stream
  
  repeat {
    fq=yield(stream)
    if (length(fq)==0)
      break
    ## get fastq chunk
    
    ###TRIM
    fq = narrow(fq,start=5,end=70)
    ## trim the first 4 and the last 2 bases

    ####FILTER
    fq = fq[myFilter(fq)]
    ## remove reads that:
      ## contain more than 1 N
      ## have median quality < 30
    
    writeFastq(fq, destination, mode="a")
    ## Append to fastq file
  }
}

## ----ShortRead_ApplyFilter,eval=F,echo=TRUE,cache=TRUE-------------------
## FilterAndTrim(fqFiles[1],
##               destination=file.path(getwd(),"FilteredFastq.fastq"))
## FilteredFastq=readFastq("FilteredFastq.fastq")
## FilteredFastq

## ----SAMBAM_title, eval=FALSE,echo=FALSE, message=FALSE------------------
## ##---------------------------------------------------------##
## ## Working with SAM/BAM files
## ##---------------------------------------------------------##

## ----Rsamtools_lib,eval=FALSE,echo=FALSE---------------------------------
## require(Rsamtools, quietly = T)

## ----GenomicAlignments_lib,eval=FALSE,echo=FALSE-------------------------
## require(GenomicAlignments, quietly = T)

## ----ToolsForSAMBAM_title, eval=FALSE,echo=FALSE, message=FALSE----------
## ##----------##
## ## Tools for SAM/BAM files
## ##--------- ##

## ----Importing_BAM_title, eval=FALSE,echo=FALSE, message=FALSE-----------
## ##----------##
## ## Importing BAM files
## ##--------- ##

## ----sr_bamfile_path,cache=TRUE------------------------------------------
sr_bamFile=untreated1_chr4() # from passilaBamSubset package 

## ----build_BAM_index,eval=FALSE,echo=TRUE--------------------------------
## indexBam(sr_bamFile)

## ----ScamBamParam_define,eval=TRUE,echo=TRUE,cache=TRUE------------------
which = RangesList("chr2L"=IRanges(7000,10000),
                   "chr4"=IRanges(c(75000,1190000),c(85000,1203000)))

scanBamWhat() #available fields
what = c("rname","strand","pos","qwidth","seq")
flag=scanBamFlag(isDuplicate=FALSE)
param=ScanBamParam(which=which,what=what,flag=flag)

## ----srBAM_import_scanBam,cache=TRUE-------------------------------------
mysrbam=scanBam(sr_bamFile,param=param)
class(mysrbam)
names(mysrbam)
sapply(mysrbam,sapply,length)["rname",] #number of imported reads

## ----srBAM_import_readGAlignments,cache=TRUE-----------------------------
mysrbam2=readGAlignments(sr_bamFile,
                         param=ScanBamParam(which=which,
                                            what="seq",
                                            flag=flag))
mysrbam2[1:2]

## ----srBAM_import_readGAlignments_v2,cache=TRUE--------------------------
mysrbam2=readGAlignments(sr_bamFile,
                         param=ScanBamParam(which=which))
mysrbam2[1:2]

## ----GAlignments_accessors,cache=TRUE------------------------------------
head(start(mysrbam2))
head(width(mysrbam2))
seqnames(mysrbam2)
cigar(mysrbam2)[1:3]
head(njunc(mysrbam2))

## ----GAlignments_to_GRanges,cache=TRUE-----------------------------------
granges(mysrbam2)[1:2]

## ----GAlignments_getJunctions,cache=TRUE---------------------------------
grglist(mysrbam2)[[1]] #only first read shown element here
junctions(mysrbam2)[[1]] #and the corresponding junctions

## ----pr_bamfile_path,cache=TRUE------------------------------------------
pr_bamFile=untreated3_chr4() # from passilaBamSubset package

## ----readGAlignmentsPairs,cache=TRUE-------------------------------------
myprbam=readGAlignmentPairs(pr_bamFile,
                              param=ScanBamParam(which=which))
myprbam[1:2]

## ----readGAlignmentsPairs_access,cache=TRUE------------------------------
myprbam[1] #first record
first(myprbam[1]) #first sequenced fragment
last(myprbam[1]) #last sequenced fragment

## ----readGAlignmentsList,cache=TRUE--------------------------------------
myprbam_list=readGAlignmentsList(pr_bamFile,
                                 param=ScanBamParam(which=which))
myprbam_list[1:2]
table(elementNROWS(myprbam_list))
summary(mcols(myprbam_list)$mate_status) #mate status as metadata

## ----Looping_on_BAM,cache=TRUE-------------------------------------------
bf = BamFile(sr_bamFile,yieldSize=100000) #create reference
open(bf) #open connection
cvg = NULL #initialize

repeat {
  chunk <- readGAlignments(bf) #loop on the BAM file
  if (length(chunk) == 0L)
    break
  chunk_cvg <- coverage(chunk)
  if (is.null(cvg)) {
    cvg <- chunk_cvg
  } else {
    cvg <- cvg + chunk_cvg
  }
}
close(bf)

cvg$chr4

## ----QAQC_on_BAM_title, eval=FALSE,echo=FALSE, message=FALSE-------------
## ##----------##
## ## Some QA/QC on aligned reads
## ##--------- ##

## ----quickBamFlagSummary,cache=TRUE--------------------------------------
quickBamFlagSummary(pr_bamFile)

## ----Coverage_title, eval=FALSE,echo=FALSE, message=FALSE----------------
## ##----------##
## ## Computing a coverage
## ##--------- ##

## ----Coverage_calc,cache=TRUE--------------------------------------------
cvg_sr=coverage(mysrbam2)
cvg_sr$chr4

## ----coverage_extract_Views,cache=TRUE-----------------------------------
#convert imported intervals to GRanges:
which_chr4_gr=GRanges(seqnames="chr4",
                      ranges=which$chr4,strand="*") 
#Get the exons
ex_chr4=exonsByOverlaps(TxDb.Dmelanogaster.UCSC.dm3.ensGene,which_chr4_gr) 
head(Views(cvg_sr$chr4,ranges(ex_chr4))) #extract Views

## ----coverage_By_strand,cache=TRUE---------------------------------------
coverage(mysrbam2[strand(mysrbam2)=="-"])$chr4

## ----Get_ChIPseq_data,cache=TRUE-----------------------------------------
# library(MMDiffBamSubset)
ChIPex_path = WT.AB2() #from MMDiffBamSubset package
ChIP_ga=readGAlignments(ChIPex_path,
                        param=ScanBamParam(
                          which=GRanges(seqnames="chr1",
                                        ranges=IRanges(3e6,5e6),
                                        strand="*")))

## ----coverage_shift_reduce,cache=TRUE------------------------------------
coverage(ChIP_ga)$chr1
#directly in coverage (!shift is not strand-aware):
coverage(ChIP_ga,shift=150)$chr1 
#now "strand-aware" shift:
coverage(ChIP_ga,shift=150*as.numeric(paste(strand(ChIP_ga),"1",sep="")))$chr1 
#resize via a GRanges (strand-aware):
coverage(resize(granges(ChIP_ga),300))$chr1 

## ----Finding_peaks_title, eval=FALSE,echo=FALSE, message=FALSE-----------
## ##----------##
## ## Finding peaks in read coverage
## ##--------- ##

## ----Slice_on_coverage,cache=TRUE----------------------------------------
cvg_H3K4me3=coverage(resize(granges(ChIP_ga),300))
slice(cvg_H3K4me3,lower=20)$chr1

## ----CountingReads_title, eval=FALSE,echo=FALSE, message=FALSE-----------
## ##----------##
## ## Counting reads / read summarization
## ##--------- ##

## ----GetExonsByGene_chr4,cache=TRUE--------------------------------------
exbygn_chr4=subsetByOverlaps(exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                                     by="gene"),ex_chr4)

## ----summarizeOverlaps_ex1,cache=TRUE------------------------------------
count_res=summarizeOverlaps(exbygn_chr4, mysrbam2, mode="Union")
count_res
assays(count_res)$counts

## ----summarizeOverlaps_ex2,cache=TRUE------------------------------------
count_res2=summarizeOverlaps(exbygn_chr4,
                             c(sr_bamFile,pr_bamFile),
                             mode="Union")

assays(count_res2)$counts

## ----summarizeOverlaps_ex3,cache=TRUE------------------------------------
assays(summarizeOverlaps(exbygn_chr4,myprbam,mode="Union"))$counts
assays(summarizeOverlaps(exbygn_chr4,myprbam,mode="Union", #ignore strand
                         ignore.strand=T))$counts

assays(summarizeOverlaps(exbygn_chr4,mysrbam2,mode="Union"))$counts
#resize the reads to 100:
assays(
    summarizeOverlaps(exbygn_chr4,mysrbam2,mode="Union", 
                      preprocess.reads=function(x){
                          resize(granges(x),100)
                          })
    )$counts

## ----Annotation_title, eval=FALSE,echo=FALSE, message=FALSE--------------
## ##---------------------------------------------------------##
## ## More annotation packages in Bioconductor
## ##---------------------------------------------------------##

## ----AnnotationDbi_lib,eval=FALSE,echo=FALSE-----------------------------
## require(AnnotationDbi)

## ----org.Dm.eg.db_lib,eval=FALSE,echo=FALSE------------------------------
## require(org.Dm.eg.db)

## ----drosophila2_ChipDb_lib,eval=FALSE,echo=FALSE------------------------
## require(drosophila2.db)
## require(drosophila2probe)
## require(drosophila2cdf)

## ----hom.Dm.inp.db_lib,eval=FALSE,echo=FALSE-----------------------------
## require(hom.Dm.inp.db)

## ----GO.db_lib,eval=FALSE,echo=FALSE-------------------------------------
## require(GO.db)

## ----TxDb.Dmelanogaster.UCSC.dm3.ensGene_lib2,eval=FALSE,echo=FALSE------
## require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

## ----TypesOfAnnotPackage_title, eval=FALSE,echo=FALSE, message=FALSE-----
## ##----------##
## ## Types of annotation packages
## ##--------- ##

## ----AccessingAnnotations_title, eval=FALSE,echo=FALSE, message=FALSE----
## ##----------##
## ## Accessing annotations
## ##--------- ##

## ----Annotation_OrgDb_access,eval=-c(1,3,4),echo=TRUE,cache=TRUE---------
library(org.Dm.eg.db)
columns(org.Dm.eg.db)
help("PATH")
keytypes(org.Dm.eg.db) #same as columns(org.Dm.eg.db) in this case
uniKeys = keys(org.Dm.eg.db,keytype="UNIPROT")[c(5,6,24)]

cols = c("SYMBOL","GO")
select(org.Dm.eg.db,
       keys=uniKeys[1:2],
       columns=cols,
       keytype="UNIPROT")

## ----Annotation_GOdb,eval=-1,echo=TRUE,cache=TRUE------------------------
library(GO.db)
mygos=c("GO:0000791","GO:0002165", "GO:0008270")
select(GO.db,
       columns=columns(GO.db)[2:4],
       keys=mygos,
       keytype="GOID")

## ----Annotation_GOdb_toTable,eval=-1,echo=TRUE,cache=TRUE----------------
ls("package:GO.db")
toTable(GOTERM)[1:3,2:4]

## ----Annotation_OrgDb_keySearch,cache=TRUE-------------------------------
keys(org.Dm.eg.db, 
     keytype="SYMBOL", 
     pattern="EcR")

select(org.Dm.eg.db,
       keys=c("EcR","DopEcR"),
       columns=c("ENTREZID","FLYBASE","GENENAME"),
       keytype="SYMBOL")

## ----Annotation_ChIPDb,eval=-c(1,3),echo=TRUE,cache=TRUE-----------------
library(drosophila2.db)
ls("package:drosophila2.db")[1:8]
drosophila2.db #provides information on the underlying database
columns(drosophila2.db)[1:7]
select(drosophila2.db,columns=c("ENTREZID","SYMBOL","ENSEMBL"),
       keys=c("1639797_at","1627097_at","1628020_at"),keytype="PROBEID")

## ----Annotation_probeAndcdf,eval=-c(1,2),echo=TRUE,cache=TRUE------------
library(drosophila2probe)
library(drosophila2cdf)
drosophila2probe[1:2,] #a data frame with probe sequences
ls("package:drosophila2cdf")

## ----Annotation_InparanoidDb,eval=-1,echo=TRUE,cache=TRUE----------------
library(hom.Dm.inp.db)
select(hom.Dm.inp.db,
       columns=c("HOMO_SAPIENS","CULEX_PIPIENS"),
       keys=c("FBpp0084497", "FBpp0077213"),
       keytype="DROSOPHILA_MELANOGASTER")
select(hom.Dm.inp.db,
       columns=c("HOMO_SAPIENS","DROSOPHILA_MELANOGASTER"),
       keys=c("CPIJ014347","CPIJ005780"),
       keytype="CULEX_PIPIENS")

## ----Annotation_TxDb,cache=TRUE------------------------------------------
select(TxDb.Dmelanogaster.UCSC.dm3.ensGene,
       columns=c("TXID","TXCHROM","TXSTRAND","TXSTART","TXEND"),
       keys=c("FBgn0039183","FBgn0264342","FBgn0030583"),
       keytype='GENEID')

## ----ImportExport_title, eval=FALSE,echo=FALSE, message=FALSE------------
## ##---------------------------------------------------------##
## ## Import/export of genomic data
## ##---------------------------------------------------------##

## ----rtracklayer_lib,eval=FALSE,echo=FALSE-------------------------------
## require(rtracklayer, quietly = T)

## ----Gviz_lib,eval=FALSE,echo=FALSE--------------------------------------
## require(Gviz, quietly = T)

## ----AnnotationHub_lib,eval=FALSE,echo=FALSE-----------------------------
## require(AnnotationHub, quietly = T)

## ----biomaRt_lib,eval=FALSE,echo=FALSE-----------------------------------
## require(biomaRt, quietly = T)

## ----SRADb_lib,eval=FALSE,echo=FALSE-------------------------------------
## require(SRAdb, quietly = T)

## ----GEOquery_lib,eval=FALSE,echo=FALSE----------------------------------
## require(GEOquery, quietly = T)

## ----rtracklayer_title, eval=FALSE,echo=FALSE, message=FALSE-------------
## ##----------##
## ## The rtracklayer package
## ##--------- ##

## ----Path_to_ExampleFiles,cache=TRUE-------------------------------------
bamExFile_path=system.file(package="Gviz","extdata","test.bam")
gff3ExFile_path=system.file(package="Gviz","extdata","test.gff3")
gtfExFile_path=system.file(package="Gviz","extdata","test.gff2")
bedExFile_path=system.file(package="Gviz","extdata","test.bed")
wigExFile_path=system.file(package="Gviz","extdata","test.wig")
bedGraphExFile_path=system.file(package="Gviz","extdata","test.bedGraph")

## ----rtracklayer_import--------------------------------------------------
head(import(bedExFile_path))
head(import(wigExFile_path)) #binned at 300bp
Rle(rep(import(wigExFile_path)$score,each=300)) #convert to Rle

## ----AnnotationHub_title, eval=FALSE,echo=FALSE, message=FALSE-----------
## ##----------##
## ## The AnnotationHub package
## ##--------- ##

## ----AnnotationHub_lib2,eval=FALSE,echo=TRUE-----------------------------
## library(AnnotationHub)

## ----AnnotationHub_startSession,eval=TRUE,echo=TRUE,cache=TRUE-----------
ah=AnnotationHub()

## ----AnnotationHub_Explore_ah,eval=TRUE,echo=TRUE,cache=TRUE-------------
annot_ah=mcols(ah) #Informations on the different records
table(annot_ah$rdataclass) #type of files that can be retrieved
table(annot_ah$dataprovider)[1:6] #providers

## ----AnnotationHub_orgDb_example,eval=TRUE,echo=TRUE,cache=TRUE----------
query(ah,"orgDb") #search for available orgDb packages
query(ah,c("orgDb","Arabidopsis")) #those for Arabidopsis only
newAt=ah[["AH57965"]] #retrieve the orgDb object
keytypes(newAt) #explore the object
select(newAt,
       keys="AT1G01010",
       keytype="TAIR",
       columns=c("REFSEQ","ENTREZID","GO"))

## ----AnnotationHub_EpigenomeRoadmap_example,eval=TRUE,echo=TRUE,cache=TRUE----
epiFiles=query(ah,"EpigenomeRoadMap")
epiFiles
unique(epiFiles$species) # sanity check
unique(epiFiles$genome) # sanity check
table(epiFiles$sourcetype) #types of files available
#more precise description of the files available:
head(sort(table(epiFiles$description), decreasing=TRUE)) 
#a more precise query:
query(ah , c("EpigenomeRoadMap","H3K36ME3","broadPeak","liver")) 
k36Peaks=ah[["AH29351"]] #retrieve the data
k36Peaks # a GRanges object
metadata(k36Peaks)

## ----AnnotationHub_LiftOver_example,eval=TRUE,echo=TRUE,cache=TRUE-------
query(ah,c("dm3","dm6","chainfile")) #search for a chain file
chain=ah[["AH15105"]] #retrieve the chain file
genes #Drosophila genes
liftOver(genes,chain) #new coordinates

## ----AnnotationHub_versions,eval=TRUE,echo=TRUE,cache=TRUE---------------
possibleDates(ah)[1:10] # available dates for the Hub
snapshotDate(ah) #date currently in use (can be changed using <-)

## ----AnnotationHub_cacheLocation-----------------------------------------
hubCache(ah)
hubUrl(ah)

## ----AnnotationHub_FDisplay_ah,eval=FALSE,echo=TRUE----------------------
## d=display(ah)

## ----biomaRt_title, eval=FALSE,echo=FALSE, message=FALSE-----------------
## ##----------##
## ## The biomaRt package
## ##--------- ##

## ----biomaRt_listMarts,eval=-1,cache=TRUE--------------------------------
library(biomaRt)
listMarts(host="www.ensembl.org")

## ----biomaRt_selectMart,cache=TRUE---------------------------------------
ens=useMart('ENSEMBL_MART_ENSEMBL',host='www.ensembl.org')
ens
listDatasets(ens)[1:5,]
rattus=useMart('ENSEMBL_MART_ENSEMBL',
               host='www.ensembl.org',
               dataset='rnorvegicus_gene_ensembl')

## ----biomaRt_GOI,cache=TRUE----------------------------------------------
head(keytypes(rattus)) #see also listFilters(rattus)
head(columns(rattus)) #see also listAttributes(rattus)

## ----biomaRt_keysFUN,cache=TRUE------------------------------------------
keys(rattus,keytype="chromosome_name")[1:5]

## ----biomaRt_getBM,cache=TRUE--------------------------------------------
goi=c('ENSRNOG00000012586','ENSRNOG00000018113') #genes of interest
getBM(attributes=c('ensembl_gene_id', 'strand',
                   'chromosome_name','start_position','end_position'), 
                    filters = 'ensembl_gene_id', 
                    values = goi, mart = rattus)

## ----biomaRt_select,cache=TRUE-------------------------------------------
select(rattus,keys=goi,keytype='ensembl_gene_id',
       columns=c('ensembl_gene_id', 'strand',
                 'chromosome_name','start_position','end_position'))

## ----GEOquery_SRAdb_title, eval=FALSE,echo=FALSE, message=FALSE----------
## ##----------##
## ## The GEOquery and SRAdb packages
## ##--------- ##

## ----GEOquery_rawdataDownload,eval=FALSE,echo=TRUE-----------------------
## library(GEOquery)
## RawGSE13149=getGEOSuppFiles('GSE13149') #!! large files !!

## ----GEOquery_getFiles,eval=FALSE,echo=TRUE------------------------------
## gse13149=getGEO('GSE13149')
## show(gse13149) #Here only one Expression Set
## GEOeset=gse13149[[1]] #get the Expression Set

## ----SRAdb_getSRAdbFile,eval=FALSE---------------------------------------
## library(SRAdb)
## sqlfile = 'SRAmetadb.sqlite'
## if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile() #large file!!

## ----SRAdb_connect,eval=FALSE--------------------------------------------
## sra_con = dbConnect(SQLite(),sqlfile)

## ----SRADb_exploreContent,eval=FALSE-------------------------------------
## dbListTables(sra_con) #Tables available in the database
## dbListFields(sra_con,"study") #Fields for the study table
## colDesc=colDescriptions(sra_con=sra_con) #Description of the fields

## ----SRADb_Query,eval=FALSE----------------------------------------------
## rs = dbGetQuery(sra_con,"select * from study limit 3")
## rs[, 1:3]
## ##   study_ID              study_alias study_accession
## ## 1        1                DRP000001       DRP000001
## ## 2        2                DRP000002       DRP000002
## ## 3        3 DLD1_normoxia_nucleosome       DRP000003

## ----SRADb_sraConvert,eval=FALSE-----------------------------------------
## sraConvert( 'SRP001007', sra_con = sra_con )
## ##       study submission    sample experiment       run
## ## 1 SRP001007  SRA009276 SRS004650  SRX007396 SRR020740
## ## 2 SRP001007  SRA009276 SRS004650  SRX007396 SRR020739

## ----SRADb_getSRA,eval=FALSE---------------------------------------------
## rs = getSRA( search_terms = "RNASeq",
##               out_types = c('study'), sra_con )
## dim(rs)
## ## [1] 1172   12
## head(rs[,1:2])
## ##   study_alias     study
## ## 1   DRP000366 DRP000366
## ## 2   PRJDB2653 DRP000375
## ## 3   PRJDB2395 DRP000376
## ## 4   PRJDB2135 DRP000535
## ## 5   PRJDB2122 DRP000536
## ## 6   PRJDB2739 DRP000537

## ----SRADb_listSRAfile,eval=FALSE----------------------------------------
## listSRAfile( c("SRX000122"), sra_con, fileType = 'sra' )
## #or: getSRAinfo("SRX000122",sra_con,sraType="sra")

## ----SRAdb_disconnect,eval=FALSE-----------------------------------------
## dbDisconnect(sra_con)

## ----Visualization_title, eval=FALSE,echo=FALSE, message=FALSE-----------
## ##---------------------------------------------------------##
## ## Visualization of genomic data
## ##---------------------------------------------------------##

## ----Gviz_lib2,eval=FALSE,echo=FALSE-------------------------------------
## require(Gviz, quietly = T)

## ----Intro_visualization_title, eval=FALSE,echo=FALSE, message=FALSE-----
## ##----------##
## ## Introduction to R graphics for genomic data
## ##--------- ##

## ----Gviz package_title, eval=FALSE,echo=FALSE, message=FALSE------------
## ##----------##
## ## The Gviz package
## ##--------- ##

## ----Gviz_ROI,cache=TRUE-------------------------------------------------
ROI=GRanges(seqnames="chr4",ranges=which$chr4[1],strand="*")

## ----Gviz_TATAAAgr,cache=TRUE--------------------------------------------
TATAAA_on_ROI=shift(union(ranges(matchPattern('TATAAA',
                               subseq(Dmelanogaster$chr4,
                                      start=start(ROI),
                                      end=end(ROI)))),
                 ranges(matchPattern('TTTATA',
                              subseq(Dmelanogaster$chr4,
                                     start=start(ROI),
                                     end=end(ROI))))),
                 start(ROI))
TATAAAgr=GRanges(seqnames="chr4",
                        ranges=TATAAA_on_ROI,
                        strand="*")

## ----Gviz_AnnotationTrack,cache=TRUE-------------------------------------
atrack=AnnotationTrack(TATAAAgr,name="TATAAA motif")

## ----Gviz_GenomeAxisTrack,cache=TRUE-------------------------------------
gtrack=GenomeAxisTrack()

## ----Gviz_IdeogramTrack,eval=FALSE,cache=TRUE----------------------------
## itrack=IdeogramTrack(genome="dm3", chromosome="chr4")
## #I'm having some unresolved issues with these ideograms
## #so I don't use them below

## ----Gviz_gb1,fig.width=10, fig.height=4,out.width="0.7\\textwidth",fig.pos="!h",fig.align="center",fig.cap='Genome axis and TATAAA Annotationtracks.',cache=TRUE----
plotTracks(list(gtrack,atrack)) #add itrack if possible

## ----Gviz_GeneRegionTrack,cache=TRUE-------------------------------------
grtrack=GeneRegionTrack(TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                        start=start(ROI),
                        end=end(ROI),
                        genome="dm3",chromosome="chr4",
                        name="Gene Model")

## ----Gviz_SequenceTrack,cache=TRUE---------------------------------------
strack=SequenceTrack(Dmelanogaster, chromosome="chr4")

## ----Gviz_gb2,fig.width=10, fig.height=4,out.width="0.74\\textwidth",fig.pos="!h",fig.align="center",fig.cap='Visualizing Gene and Sequence Tracks',cache=TRUE----
plotTracks(list(gtrack,atrack,grtrack,strack))

## ----Gviz_gb3_zooms,fig.width=10, fig.height=4,fig.align="center",fig.pos="!h",out.width="0.74\\textwidth",fig.show='asis',fig.cap='Zomming in and out',cache=TRUE----
plotTracks(list(gtrack,atrack,grtrack),
           extend.left = 0.5, extend.right = 10000)

## ----Gviz_gb3_zooms2,fig.width=10, fig.height=4,fig.align="center",fig.pos="!h",out.width="0.74\\textwidth",fig.show='asis',fig.cap='Zomming in',cache=TRUE----
plotTracks(list(gtrack,atrack,strack),
           from = 79900, to = 80100)

## ----Gviz_gb4_zoom,fig.width=8, fig.height=4,out.width="0.74\\textwidth",fig.align="center",fig.pos="!h",fig.cap='Zomming in some more to read the sequence',cache=TRUE----
plotTracks(list(gtrack,atrack,strack),
           from = 79975, to = 80015)

## ----Gviz_AlignmentsTrack,cache=TRUE-------------------------------------
altrack=AlignmentsTrack(pr_bamFile, isPaired=TRUE)

## ----Gviz_gb5_reads,fig.width=10,fig.height=5,fig.align="center",fig.pos="!h",out.width="0.72\\textwidth",fig.show='asis',fig.cap='Visualizing Alignment Tracks.',cache=TRUE----
plotTracks(list(gtrack,atrack,
                grtrack,altrack)) #use type="coverage" to see only the coverage

## ----Gviz_gb5_reads2,fig.width=10,fig.height=5,fig.align="center",fig.pos="!h",out.width="0.72\\textwidth",fig.show='asis',fig.cap='Visualizing Alignment Tracks.',cache=TRUE----
plotTracks(list(gtrack,atrack,grtrack,altrack,strack),
           from=72000,to=73500,
           col.mates="purple",
           col.gaps="orange")

## ----Gviz_DataTrack,cache=TRUE-------------------------------------------
set.seed(255)
lim <- c(start(ROI), end(ROI))
coords <- sort(c(lim[1], sample(seq(from=lim[1], to=lim[2]), 99), lim[2]))
c1=runif(100, min=-10, max=8)
c2=runif(100, min=-10, max=8)
dat=GRanges(seqnames="chr4",strand="*",
            ranges=IRanges(start=coords[-length(coords)], end=coords[-1]),
            ctrol1=c1,ctrol2=c2,
            treated1=c1+rnorm(100,3),
            treated2=c2+rnorm(100,3))
genome(dat)="dm3"
dtrack <- DataTrack(dat,name="Uniform")

## ----Gviz_gb6_DataTrack,fig.width=8, fig.height=5,fig.align="center",out.width="0.7\\textwidth",fig.pos="!h",fig.show='hold',fig.cap='Visualizing Data Tracks.',cache=TRUE----
plotTracks(list(gtrack, atrack, grtrack, dtrack),
           from=lim[1], to=lim[2], type=c("a","p"),
           groups=rep(c("ctrol","treated"),each=2))

## ----SessionInfo,eval=TRUE,echo=FALSE,results="asis"---------------------
toLatex(sessionInfo())

