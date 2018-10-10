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
## Ex02: 
##-------------------------------

