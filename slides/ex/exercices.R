##-------------------------------
## Ex01: Practice with your own DNA sequence
##-------------------------------

library(Biostrings)

set.seed(123)
myseq <- paste(sample(DNA_ALPHABET[1:4],
                      size = 50,
                      replace = TRUE),
               collapse = "")

myseq <- DNAString(myseq)

reverseComplement(myseq)

100 * alphabetFrequency(myseq, baseOnly=TRUE, as.prob = TRUE)[1:4]

letterFrequency(myseq, "CG", as.prob = TRUE)

myseq[11:20]
subseq(myseq, 11, 20)

toString(myseq[11:20])
as.character(myseq[11:20])

Views(myseq, start=seq(1,50, by=10), width=5)


##-------------------------------
## Ex02: Motif pipeline example
##-------------------------------

#Get a GRanges of peaks via AnnotationHub:
library(AnnotationHub)
library(rtracklayer)
library(Biostrings)
ah <- AnnotationHub()

#get some files with peaks for CTCF in human HepG2 cells
ctcfFiles <- query(ah, c("ENCODE", "hg19", "CTCF", "Hepg2", "narrowPeak"))
# Available metadata
colnames(mcols(ctcfFiles))

AllctcfPeaks <- lapply(names(ctcfFiles),
                       function(nn) ah[[nn]])
AllctcfPeaks <- as(AllctcfPeaks,
                   "GRangesList")

#For each peak count the number of datasets in which a peak overlapping this peak appears
NumOVL <- lapply(AllctcfPeaks,
                 function(x) {countOverlaps(x, AllctcfPeaks)})
#Number of peaks appearing in all 7 datasets:
sapply(NumOVL, function(x){sum(x==7)})

# Select the peaks "present" in all datasets
SelectedPeaks <- unlist(
                     as(mapply(function(peaks, ovl) {peaks[ovl==7]},
                               AllctcfPeaks,
                               NumOVL),
                        "GRangesList"))

# Reduce (we get ~34K peaks)
SelectedPeaks <- reduce(SelectedPeaks)

# To limit the computation time, we select only the peaks on the autosomes that are <250bp and resize them to 200bp
## Note that it would be smarter to filter on the significance of the peaks (available in the initial GRanges)
SelectedPeaks <- resize(SelectedPeaks[width(SelectedPeaks)<250],
                        width=200, fix="center")
SelectedPeaks <- SelectedPeaks[seqnames(SelectedPeaks) %in% paste0("chr",1:22)]
#We get ~4500 peaks

#Get the corresponding sequences
library(BSgenome.Hsapiens.UCSC.hg19)
SelSeq <- getSeq(Hsapiens, SelectedPeaks)

# Use rGADEM to search for motifs (note that GADEM could use the score of ChIP-seq peaks: see ?GADEM)
library(rGADEM)
gadem <- GADEM(SelSeq, verbose=1, genome=Hsapiens) #using an unseeded analysis (!! very long !!)
## The object is saved is /ex folder.
## Import using: gadem <- readRDS("/ex/gadem.rds")

#Take a look at the 6 motif founds
consensus(gadem) #consensus motifs
nOccurrences(gademObj) #Number of occurences of the 6 motifs in the sequences
#get the PFMs:
gadem_PFMs <- lapply(1:nMotifs(gadem),
                     function(num){gadem@motifList[[num]]@pwm})
names(gadem_PFMs) <- names(gadem)
#convert them to motifStack pfm and plot:
library(motifStack)
gadem_pfm <- lapply(names(gadem_PFMs),
                    function(x) {new("pfm",
                                     mat=gadem_PFMs[[x]],
                                     name = x)})
motifStack(gadem_pfm)


#We could use the GADEM seeded analysis in which we input a motif that is likely to be found but the computation is still very long:
library(TFBSTools)
library(JASPAR2018)
#~ CTCF <- TFBSTools::getMatrixByID(JASPAR2018, "MA0531")
#~ ctcf.pfm <- apply(Matrix(CTCF), 2 , as.numeric)
#~ rownames(ctcf.pfm) <- rownames(Matrix(CTCF))
#~ gademSeeded <- GADEM(SelSeq, verbose = 1, genome = Hsapiens, Spwm = list(ctcf.pfm))
#~ saveRDS(gademSeeded, "/DATA/WORK/Enseignement/BioC_For_NGS/slides/ex/gademSeeded.rds")

# Use MEME via TFBSTools (requires MEME to be installed)
library(TFBSTools)
#~ memeMotif <- runMEME(SelSeq, binary="meme", arguments=list("-nmotifs"=3))



##-------------------------------
## Ex03: Compare indexes
##-------------------------------

library(Biostrings)

indx <- DNAStringSet( c("ATCACG", "CGATGT", "TTAGGC", "TGACCA",
                        "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA",
                        "GATCAG", "TAGCTT", "GGCTAC", "CTTGTA",
                        "CGGCTA", "TCCGCG", "ATGTCA", "AGCGAT"))

## Get the edit distances between the indexes:
indxDist <- stringDist(indx)
indxDist <- as.matrix(indxDist)
diag(indxDist) <- NA
apply(indxDist, 1, min, na.rm=TRUE)
#We can remove the indexes with a minimum edit distance of 2

## Note that the stringDist with method = "levenshtein" also considers indels which might be unlikely
## If we want to only count the number of substitutions only, we can use:
subst <- as.matrix(
         stringDist(indx, 
                    method = "substitutionMatrix", 
                    substitutionMatrix = 1-nucleotideSubstitutionMatrix()[1:4,1:4]))
diag(subst)=NA
#Get the minimum number of substitution between indexes
apply(subst, 1, min, na.rm=T)
#Here we see that all indexes have at least 3 substitutions between them

## Select the indexes
selindx <- indx[apply(indxDist, 1, min, na.rm=TRUE)>2]

## Get the consensusMatrix
consmat <- consensusMatrix(indx)[1:4,]
#Normalize the matrix in frequencies:
consmatnorm <- t(t(consmat)/colSums(consmat))
#image the proportion at each position:
image(t(consmatnorm[4:1,]),col=colorRampPalette(c("darkgreen","white","darkred"))(11),
      axes=F, zlim=c(0,0.5))
axis(side=2,at=c(0,0.33,0.66,1), label=rev(c("A","C","G","T")), tick=F, las=2)
axis(side=3,at=seq(0,1,length.out=ncol(consmat)), label=1:ncol(consmat), tick=F)


##-------------------------------
## Ex03: 
##-------------------------------



