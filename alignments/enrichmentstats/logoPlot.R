require(Biostrings)
require(ggplot2)
require(ggseqlogo)

svg(filename="2xb7_inactives.svg", 
    width=5, 
    height=4, 
    pointsize=12)
fastaFile="../trees/2xb7_manning_blossum90_manning.fastainactives.fasta"
#read.fasta(file=fastaFile)
fa=readAAStringSet(filepath=fastaFile)
#fa2=scan(fastaFile, character(), quote = "")
#fa2[!grepl('^>', fa2)]

ggseqlogo(data.frame(fa)$fa )
dev.off()