source("~/local/DEAFunctions.R")

mkMatx4enh<-function(dir){

   #Get the files-REMEMBER to fill in the paste arguments with the <samout> option set for htseq once you run it-DONE
   files<-list.files(path=dir, pattern= glob2rx(paste("*",".counts", sep="")), full.names=TRUE)
   files<-files[grepl("*.(plus|minus).counts", files)]
   #Save the names (cells and conditions) for the files-REMEMBER to fill it in-DONE
   #sub("pattern", "replacement", x)
   names<-gsub("\\..*enhancer", "", files)
   names<-gsub(".*/", "", names)

   #Read the files and create a list of tables
   counts_list<-lapply(files, read.table, row.names=1)
   #Initialize a matrix
   counts_matrix<- NULL
   #Iterate through the elements of the list (i.e. the tables)

   for(i in 1:length(counts_list)){
      #for(i in 1:length(counts_list)){
      #Stop and produce an error message indicating the first of the elements of which were not true if not all tables have the same genes (rownames)
      stopifnot(all(rownames(counts_list[[i]])== rownames(counts_list[[1]])))
      #make a matrix. Our table only has two columns, with the first one (V1) being NULL because contains the genes and the V2 containing the counts for each gene
      counts_matrix<-cbind(counts_matrix, counts_list[[i]]$V2)
      }

   #Get the gene names
   rownames(counts_matrix)<-rownames(counts_list[[1]])
   colnames(counts_matrix)<-names

   #Removing the last 5 lines of the HTSeq-Count output
   counts_matrix<-subset(counts_matrix, grepl("E", rownames(counts_matrix)))
   return(counts_matrix)
}

draw.enhanc.heatmap <- function(matrix, list, title, rsc){
   library("gplots")
   library("RColorBrewer")
   require(pheatmap)

   m <- matrix[list,]
   width <- length(list) / 5
    if (length(list) > 1) {
       plot.data <- as.matrix(log(m+1))
       plot.data <- plot.data-rowMeans(plot.data)

       plot.data <- t(plot.data)

       c.ln <- 254
       c <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(c.ln)
       #cool.col <- colorRampPalette(c('azure', 'darkcyan'))(255)
       ann.col <- c('grey', 'firebrick')
       names(ann.col) <- c('N','Y')
       ann.expr<- c('red3', 'royalblue3', 'green4')
       names(ann.expr) <- c('up', 'down', 'notDE/notGenic')
       ann.l <- list(ann.col, ann.col, ann.expr)
       names(ann.l) <- colnames(rsc)

#       r <- ceiling(max(abs(min(plot.data)), max(plot.data)))
#       r <- c(seq(from=-r, to=-0.0000000001, length.out=c.ln/2), 0, seq(from=0.0000000001, to=r, length.out=c.ln/2))

       #plot
       if (nrow(plot.data) > 1) {
         pheatmap(plot.data,
                  color = c,
                  kmeans_k = NA, breaks = NA, border_color = NA,
                  cellwidth = 6, cellheight = 18, scale = "none", #could you scale= row if I don't use rowMeans before
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  clustering_distance_rows = "euclidean",
                  #clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  #treeheight_row = ifelse(cluster_rows, 21, 0),
                  #treeheight_col = ifelse(cluster_cols, 21, 0),
                  legend = TRUE, legend_breaks = NA, legend_labels = NA,
                  annotation = rsc, annotation_colors = ann.l,
                  #annotation_legend = TRUE, drop_levels = TRUE,
                  show_rownames = T, show_colnames = T, main = title,
                  fontsize = 8, #fontsize_row = fontsize,
                  #fontsize_col = fontsize, display_numbers = F,
                  number_format = "%.2f",
                  #fontsize_number = 0.8 * fontsize, filename = NA,
                  width = NA, height = NA)
       }
   }
}

myHeatmap.enh <-function(matrix, list, name, rsc){
   library("gplots")
   library("RColorBrewer")
   require("pheatmap")

   list <- list[list %in% rownames(matrix)]
   width <- length(list) / 5
    if (length(list) > 1) {
       list.matrix <- matrix[list,]
       plot.data <- as.matrix(log(list.matrix+1))
       plot.data <- plot.data-rowMeans(plot.data)

       plot.data <- t(plot.data)

       # colours
       heatcol <- colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(255)
       #plot
       if (nrow(plot.data) > 1) {
          pdf(paste("Heatmap_", name,".pdf", sep= ""), paper="a4", height= max(width,5), width= 8)
          intergenic= unique(rsc[,"intergenic"])
          names(intergenic) <- unique(rsc[,"intergenic"])
          ann_colors=list(intergenic= intergenic)
          pheatmap(plot.data2, #trace= "none", density.info= "none",
                   color=heatcol, cellwidth=10, cellheight= 7, cluster_cols=TRUE, cluster_rows=TRUE,
                   annotation_legend=TRUE, border_color=NA, clustering_distance_cols=dcols,  annotation_colors=ann_colors,
                   annotation= rsc.test, fontsize=10, fontsize_row=7, show.rownames= TRUE, show.colnames=TRUE)
          dev.off()
       }
    }
}



mk.rsc.matx <- function(m, vector){
   rsc <-NULL
   for(i in 1:length(m)){
      #vector <- gsub(unname(m[i]), names(m)[i], vector)
      rsc[grep(unname(m[i]), vector, value=FALSE)] <- names(m)[i]
   }
   return(rsc)
}


