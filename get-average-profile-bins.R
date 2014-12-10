#double check
#first.exon <- lapply(split(my.annotation$genes, names(my.annotation$genes)), function(x){ x[1,] })
#last.exon <- lapply(split(my.annotation$genes, names(my.annotation$genes)), function(x){ x[length(x),]})
##

#gr.list input should always be first, annotation= grange with names the gene names, expression=df from deseq, window=names list (names = up and down)
make.binPlots<-function(grL, anno, expr, opt){
   require(Repitools)
   require(chipseq)
   #for colMedians
   require(miscTools)
   require(RColorBrewer)

   print("***")
   print(names(grL))

   if( ! "input" %in% tolower(names(grL))){
      stop("No grange object found for input ChIP\n")
   }


   expr <- expr[expr$padj < opt['fdr'],]
   expr <- expr[order(-expr$log2FoldChange), 'log2FoldChange', drop=F]
   if(length(anno) != nrow(expr)){
      expr <- expr[rownames(expr) %in% names(anno), , drop=F]
      anno <- anno[names(anno) %in% rownames(expr), , drop=F]
   }

   #estimate the fragment size
   #fr <- sapply(grL, estimate.mean.fraglen)
   #fr <- colMedians(fr)
   #s.width to use: take the mean
   #sw <- mean(round(fr/100)*100)

   cov <- featureScores(grL, anno, up=opt['up'], down=opt['down'], freq=opt['frq'], s.width=opt['sw'])

   input.i <- grep("^input$", cov@names, value=F, ignore.case=T)
   chip.i <- grep("^input$", cov@names, value=F, ignore.case=T, invert=T)
   #should always be ip-input
   cov@scores <- list(cov@scores[[chip.i]]- cov@scores[[input.i]]) # want to subset it by name instead of index
   names(cov) <- paste(toupper(cov@names[chip.i]), toupper(cov@names[input.i]), sep='-')

   rbin <- range(expr$log2FoldChange[!is.infinite(expr$log2FoldChange) & ! is.na(expr$log2FoldChange)])
   bin <- round((rbin[2]-rbin[1])/4)
   bin <- (round(max(abs(rbin))) +1)/2
   if(bin > 9) {bin=9}

   c <- brewer.pal(9, 'Blues')[(9-bin):9]
   lwd.v <- 2
   binPlots(cov, ordering=expr, ord.label='KO-WT log2FC', plot.type="line", n.bins=bin, cols=c, lwd=lwd.v)

}
