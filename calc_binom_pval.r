
tmp_file = commandArgs()[4]

tmp = read.table(file=tmp_file, sep="\t")

pvals = dbinom(as.double(tmp[,3]),as.double(tmp[,4]),as.double(tmp[,5]))
tmp = cbind(tmp,pvals)

write.table(tmp, file=tmp_file, quote = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)

