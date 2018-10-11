##-------------------------------
## Ex01: Working with your own DNA sequence
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
## Ex02: Compare indexes
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



