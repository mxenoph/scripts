source("~/source/Rscripts/functions.R")
x <- c('org.Mm.db')
lapply(x, suppressMessages(library), character.only=T)

dir <- "/nfs/research2/bertone/user/mxenoph/hendrich/htseq_counts/de_anal"

filesGeneU <- grep("-", list.files(path=dir, pattern= glob2rx(paste("*",".txt", sep="")), full.names=TRUE), invert=TRUE, value=TRUE)
filesGeneS <- grep("-de", list.files(path=dir, pattern= glob2rx(paste("*",".txt", sep="")), full.names=TRUE), value=TRUE)

#List of lists of  all genes in the sample and de lists, here reading the table from nbinom, tables from glm will be different
listOfUniverses <- lapply(filesGeneU, function(x){table<-read.table(x, sep="\t", head=TRUE, colClasses=c(NA, rep("NULL", 8))); list<-table[,1]; return(list)})
names(listOfUniverses) <-  gsub(".*/", "", gsub("\\..*", "", filesGeneU))

listOfSel <- lapply(filesGeneS, function(x){table<-read.table(x, sep="\t", head=TRUE, colClasses=c(NA, rep("NULL", 8))); list<-table[,1]; return(list)})
names(listOfSel) <-  gsub("-de", "",gsub(".*/", "", gsub("\\..*", "", filesGeneS)))

#List of factors that describe whether a gene in the universe is interesting/de or not
GeneLists <- lapply(names(listOfUniverses),
                  function(i){	myfactor <- factor(as.integer(listOfUniverses[[i]] %in% listOfSel[[i]]));
				names(myfactor) <- listOfUniverses[[i]]; return(myfactor) })

names(GeneLists) <- names(listOfUniverses)

GOres<-lapply(names(GeneLists), function(i){GOtermEnr(GeneLists[[i]], names(GeneLists[i]), "BP", "org.Mm.eg.db")})
names(GOres) <- names(GeneLists)

pdf("GOtermsWT_EpivsKO_Epi.pdf", width= 12, height= 8)
par(mar=c(2.1, 8.1, 2.1, 8.1))
par(oma=c(1.5,1,1.5,1.5))
textplot(GOres[["WT_EpivsKO_Epi"]]$table, halign="center", valign="center")
dev.off()

save(GOres, file = paste(dir, '/', 'GOres.Rdata', sep=""))
