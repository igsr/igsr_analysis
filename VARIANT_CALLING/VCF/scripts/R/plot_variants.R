require(IdeoViz)
require(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

vcf_file = args[1] #vcf file with variants
windows_size = as.integer(args[2]) #length in bp of the windows used to slice the genome
outprefix = args[3] #prefix used for output file
height = as.numeric(args[4]) #height in pixels for plot (i.e. 9800)
chromName_cex = as.numeric(args[5]) #size for chromosome names (i.e. 0.8)

outprefix=paste(outprefix,".png",sep="")

ideo <- getIdeo("hg38")

A = read.table(vcf_file, as.is = T, sep="\t", header=F)

data = GRanges(A$V1, IRanges(start = A$V2, end = A$V3) )

mcols(data)$group1 = A$V4/windows_size

png(outprefix, width = 1800, height = height,  units = "px")

plotOnIdeo(chrom = seqlevels(data),ideoTable = ideo, values_GR = data, value_cols=colnames(mcols(data)),plotType='lines',addScale = FALSE,cex.axis = 1.0,chromName_cex = chromName_cex)

dev.off()
