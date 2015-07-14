#!/usr/bin/env Rscript

# Parsing command line arguments and create output subdirectories# {{{
library = function (...) suppressMessages(base::library(...))
select = dplyr::select
library(argparse)
library(tools)
source("~/source/Rscripts/annotation-functions.R")

parser =  ArgumentParser(description="Perform differential binding analysis")
parser$add_argument('-s', '--sheet', metavar= "file", required='True', type= "character", help= "Sample sheet must have SampleID,Condition,Treatment,Replicate,bamReads,bamControl,Peaks (last 3 are the path to the respective bed files)")
parser$add_argument('-a', '--assembly', type= "character", default='mm9', help= "Give preferred assembly e.g. mm9. Default: mm9")
parser$add_argument('-k', '--keep', type= "character", default='', help= "Factors to ignore for differential binding but plot")
parser$add_argument('-d', '--de', metavar = "path", type= "character", default='', help= "Path to look for the list of de expressed genes")
parser$add_argument('-o', '--out', metavar= "path", type= "character", default= getwd(), help= "Output directory -- all subdirectories will be created here")

args = parser$parse_args()

output_path = file.path(args$out, 'DiffBind')
plot_path = file.path(output_path, 'plots')
dir.create(plot_path, recursive= TRUE)
#}}}

# Load Packages
library(DiffBind)
library(dplyr)
# For str_replace_all function
library(stringr)
library(ggplot2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
# Color scheme for divergent colors.
divergent_colors = colorRampPalette(c('#603D71', 'white', '#A4B962'))(30)
cool_cols = colorRampPalette(c('aliceblue','darkcyan'))(100)
heat_cols = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100)

# Less intrusive defaults for heatmap.2.
heatmap.2 = function (...) 
            gplots::heatmap.2(..., trace = 'none', density.info = 'none',
                              col = heat_cols)
pheatmap = function (...) 
            pheatmap::pheatmap(..., trace = 'none', density.info = 'none',
                               #color = divergent_colors, border_color = NA,
                               color = heat_cols, border_color = NA,
                               fontsize = 10, fontsize_row = 7,
                               show_colnames = TRUE)
## Functions

# Read in and format data# {{{
import_data = function(cnf_df){

    sapply(cnf_df[['Peaks']],
           function(x){
               cmd = sprintf("awk -F \"\t\" \'NR==1 {if($0 ~ /^track/){print 1}else{print 0}}\' %s", x)
               # If track=1 file contains ucsc track line
               track = pipe(cmd, open = "r")
               # close connection
               on.exit(close(track))

               flag = scan(track, quiet=TRUE)
               if(flag==1){
                   cmd = sprintf("sed -i \'1d\' %s", x)
                   scan(pipe(cmd), quiet=TRUE)
               }
           })

    # by default dba plots the correltaion heatmap. Set it to False and print on call
    peakset = dba(sampleSheet = cnf_df, bCorPlot = FALSE)
#    print(peakset)
    return(list(cnf_df = cnf_df, peakset = peakset))
}
# }}}

# doing all the work # {{{
differential_binding = function (cnf, protein){
    tmp = import_data(cnf)
    config = tmp[['cnf_df']]
    peakset = tmp[['peakset']]
    
    # Using only this peak caller data, a correlation heatmap can be generated 
    # which gives an initial clustering of the samples
    pdf(file.path(plot_path, paste0(protein, '.pdf')))
    plot(peakset)
    peakset = dba.peakset(peakset, consensus = -DBA_REPLICATE, minOverlap = 2)

    counts = dba.count(peakset, score = DBA_SCORE_TMM_MINUS_FULL,
#    counts = dba.count(peakset, score = DBA_SCORE_READS_MINUS,
                       minOverlap = 2,
                       peaks = peakset$masks$Consensus,
#                       bRemoveDuplicates = TRUE,
                       bScaleControl = TRUE, bParallel = TRUE,
                       bCorPlot = FALSE)
    plot(counts)
    dev.off()

    contrast = dba.contrast(counts, categories = DBA_CONDITION, minMembers = 2)
    
    diff_bound = dba.analyze(contrast, method = DBA_DESEQ2,
                             bSubControl= TRUE, bFullLibrarySize = TRUE,
                             bTagwise = TRUE,
                             bParallel = TRUE,
                             bCorPlot = FALSE)

    for (index in 1:length(contrast$contrasts)){
        comparison = paste(protein,
                           contrast$contrasts[[index]]$name1,
                           contrast$contrasts[[index]]$name2, sep='-')
        pdf(file.path(plot_path, paste0(comparison, '.pdf')))

        dba.plotMA(diff_bound, contrast = index, method = DBA_DESEQ2,
                   th = 0.05, bUsePval = FALSE )

        dba.plotPCA(diff_bound, contrast = index, method = DBA_DESEQ2,
                    th = 0.05, bUsePval = FALSE,
                    label = DBA_ID)

        groups_pvals = dba.plotBox(diff_bound, contrast = index, method = DBA_DESEQ2, bAll = TRUE, pvalMethod = NULL)

        dba.plotHeatmap(diff_bound, method = DBA_DESEQ2,
                        contrast = index,
                        # Threshold will be FDR and not pval
                        bUsePval = FALSE,
                        score = DBA_SCORE_TMM_MINUS_FULL,
                        th = 0.05)

        dba.plotHeatmap(diff_bound, method = DBA_DESEQ2,
                        contrast = index,
                        # Threshold will be FDR and not pval
                        bUsePval = FALSE,
                        score = DBA_SCORE_TMM_MINUS_FULL,
                        correlations=FALSE,
                        th = 0.05)

        dba.plotHeatmap(diff_bound, method = DBA_DESEQ2,
                        contrast = index, mask = diff_bound$masks$All,
                        # Threshold will be FDR and not pval
                        bUsePval = FALSE,
                        score = DBA_SCORE_TMM_MINUS_FULL,
                        correlations=FALSE,
                        th = 0.05)
        dev.off()
        
        results =  dba.report(diff_bound, method = DBA_DESEQ2,
                              contrast = index,
                              # add count data for individual samples
                              bCounts = TRUE,
                              # only include sites with an absolute Fold value greater than equal
                              # fold= 2
                              # Threshold will be FDR and not pval
                              bUsePval = FALSE,
                              th = 0.05)

        write.table(as.data.frame(results),
                    file = file.path(output_path, paste0(comparison, '.tsv')), sep="\t", quote=FALSE, row.names=FALSE)
    }
}
# }}}

# combine db regions with expression data for nearest gene, annotation = ensembl, jan2013 etc if not in_file # {{{
combine_expression = function(comparisons, annotation = "in_file"){
    library(ChIPpeakAnno)

    # anotate db sites and plot TPM for associated genes# {{{
    per_db_region = function(db_regions, genes, deseq_results){
        annotated = annotatePeakInBatch(db_regions, AnnotationData = genes)
        annotated_df = as.data.frame(annotated)
        # Filtering, keeping only db regions in or including genes or withn 5Kb from a gene TSS or TTS
        annotated_df = annotated_df %>%
                        filter(insideFeature %in% c("inside", "includeFeature") | (insideFeature %in% c("upstream", "downstream") & shortestDistance < 5000))
        colnames(annotated_df) = gsub("^feature$", "gene_id",  colnames(annotated_df))
        keep = c(colnames(annotated_df), colnames(deseq_results)[7:20])
        keep = gsub('^strand$', 'strand.x', keep)
        keep = gsub('^start_position$', 'start_position.x', keep)
        keep = gsub('^end_position$', 'end_position.x', keep)

        annotated_df = inner_join(annotated_df, deseq_results, by = "gene_id") %>% select(one_of(keep))
        # Removes genes not expressed (padj =NA) and marks the rest as significant or not
        annotated_df = annotated_df %>% filter(!is.na(padj)) %>% mutate(DE = ifelse(padj < 0.05, 'p < 0.05', 'p >= 0.05'))

        # Scale data# {{{
        plot_data = annotated_df %>% inner_join(db_regions_df, by = 'peak') %>%
                    dplyr::select(contains('Expression'),
                                  one_of(colnames(db_regions_df)[3:ncol(db_regions_df)]), gene_id, DE)
        row_annotation = data.frame(Expression = plot_data[['DE']],
                                  row.names = paste(plot_data[['peak']], plot_data[['gene_id']], sep = "-"))
        # selecting only TPM columns and scaling based on the rowMeans for the respective experiment
        # RNA seq or ChiP-seq
        # IMPORTANT: if instead mutate_each_(funx(...), ~matches(...)) I use mutate_each(funs(...), matches(...))
        # then the selected columns are not mutated but new columns holding the correct values are created.
        # Those are names varX where X is the index of the column after subseting the table
        # Once done scaling remove the rowMeans columns
        plot_data = plot_data %>% dplyr::select(-peak, -gene_id, -DE) %>%
                    log() %>% mutate(Expression_RowMean = rowMeans(.[grep("^Expression", names(.))]),
                                     Binding_RowMean = rowMeans(.[grep("^Expression", names(.), invert = T)])) %>%
                                     mutate_each_(funs(. - Expression_RowMean), vars = ~matches("Expression.*")) %>%
                                     mutate_each_(funs(. - Binding_RowMean), vars = ~matches("^[^Expression]")) %>%
                                     select(-matches('RowMean'))

        # Id rownames are not present (and dplyr removes them, so...) the annotation row
        # in the heatmap will not work
        rownames(plot_data) = rownames(row_annotation)
        # }}}

        pheatmap(plot_data,
                 clustering_distance_rows = "correlation", scale = 'none',
                 cluster_cols = FALSE, cluster_rows = TRUE,
                 annotation_row = annotate_row, annotation_legend = T,
                 show_rownames = FALSE)

        p = ggplot(annotated_df, aes(x = log2FoldChange, y = Fold, colour = DE)) + geom_point()
        p = p + theme_classic() + xlab("log2 Fold Change of Nearest Gene") + ylab("log2 Fold Change of DB regions") + ggtitle(label)
        p = p + scale_colour_manual(values = c("#AA3929", "#8E9CA3"))
        print(p)

        # Performing a GSEA on the genes associated with a db region
        gene_scores = deseq_results %>% filter(!is.na(padj)) %>% select(gene_id, padj)
        gene_scores = with(gene_scores, structure(padj, names = as.character(gene_id)))

        # Doing GSEA on de genes and de and db genes
#        ob = get_set_enrichment(gene_scores = gene_scores, label = label)
#        ob = get_set_enrichment(gene_scores = gene_scores, bound = (select(annotated_df, gene_id) %>% .[[1]]), label = label)
    }# }}}

    # annotate genes and plot TPM and binding affinity for all de genes# {{{
    per_de_gene = function(db_regions, genes, deseq_results){
        
    }# }}}

    # Match chip-seq with rna-seq# {{{
    for( x in 1:nrow(comparisons) ) {
        db_regions = as.character(comparisons[x, 'comparisons'])
        deseq_results = as.character(comparisons[x, 'expression_data'])
        
        filename = paste0(file_path_sans_ext(basename(db_regions)),
                          "_RNAseq-",
                          file_path_sans_ext(basename(deseq_results)),
                          ".pdf")
        label = paste0("binding affinity: ",
                       file_path_sans_ext(basename(db_regions)),
                       ", expression: ",
                       file_path_sans_ext(basename(deseq_results)))
        
        db_regions = read.delim(db_regions, head = T, sep = "\t")
        # Keep TMM for replicates as I will need that downstream
        # It's not easy to keep it in granges object without knowing the names
        db_regions_df = db_regions %>%
                        mutate(peak = paste(seqnames, start, end, sep = "_")) %>%
                        dplyr::select(-seqnames, -start, -end, -width, -strand, -Conc, -Fold, -p.value, -FDR)
                        
        db_regions = with(db_regions,
                          GRanges(seqnames = seqnames,
                                  IRanges(start, end),
                                  strand = strand,
                                  Fold = Fold,
                                  FDR = FDR))
        names(db_regions) = paste(seqnames(db_regions), start(db_regions), end(db_regions), sep = "_")
        
        # IMPORTANT: use read.delim instead of read.table as quote="\"" is set by default
        # if that is not set in read.table then the file is not read through he end
        deseq_results = read.delim(deseq_results, head = T, sep = "\t")

        # If de files from Meryem, gene coordinates are in deseq files. # {{{
        # No need to load get_annotation(assembly) to get gene coordinates
        # chromosome names and strands needs to be changed
        if(annotation == "in_file"){
            source("~/source/Rscripts/granges-functions.R")
            genes = deseq_results %>% select(gene_id, chromosome_name, strand, start_position, end_position, gene_name)
            genes = with(genes, 
                         GRanges(seqnames = chromosome_name,
                                 IRanges(start_position, end_position),
                                 strand = strand,
                                 gene_id = gene_id,
                                 gene_name = gene_name))
            genes = change_seqnames(genes)
            names(genes) = genes$gene_id
            
        } else {
            #ToDO: load genes from get_annotation(assembly)
        }
        # }}}
    
        pdf(file.path(plot_path, filename), paper = 'a4')
        per_db_region(db_regions, genes, deseq_results)
        dev.off()
    }# }}}

    deseq_all = lapply(comparisons %>% select(expression_data) %>% unique() %>% .[[1]] %>% as.character(),
                       function(x){
                           read.delim(x, head = T, sep = "\t")
                       })
    names(deseq_all) = comparisons %>% select(expression_data) %>% unique() %>% .[[1]] %>% as.character()

}# }}}

# main # {{{
main = function () {
#    setwd("/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013")
#    args = list()
#    args$sheet = list.files(pattern="diff-bind.config", full.names=T)

    cnf_df = read.table(args$sheet, header=T, stringsAsFactors=F)
    proteins = levels(as.factor(cnf_df[['Factor']]))
    keep = unlist(strsplit(args$keep, ':'))
    for (x in keep){
        proteins = proteins[grep(x, proteins, invert = TRUE)]
    }

    for (p in proteins) {
        print(p)
        cnf = cnf_df %>% filter(Factor %in% c(p, keep))
        differential_binding(cnf, p)
    }

    comparisons = list.files(pattern = ".tsv", path = output_path, full.name = T)
    # Keeping only the comparisons of the time points. If design is different this won't work
    # comparisons = comparisons[grep('h.-h.', comparisons)]
    comparisons = comparisons[grep('!', comparisons, invert = T)]
    comparisons = as.data.frame(comparisons)

    tmp = comparisons %>% .[[1]] %>% as.character() %>% basename() %>% file_path_sans_ext %>% strsplit("-")
    tmp = as.data.frame(do.call(rbind, tmp)) %>% select(V2, V3)
    tmp = lapply(1:nrow(tmp), function(x) sort(tmp[x,]))
    tmp = as.data.frame(do.call(rbind, tmp))
    colnames(tmp) = c('c1', 'c2')
    comparisons = cbind(comparisons, tmp)

    expression_data = list.files(pattern = "[^de].tsv", path = args$de, full.name = T)
    expression_data = as.data.frame(expression_data)

    tmp = expression_data %>% .[[1]] %>% as.character() %>%
          basename() %>% file_path_sans_ext %>% str_replace_all("_DE|_de", "") %>% strsplit("-")
    tmp = as.data.frame(do.call(rbind, lapply(tmp, sort)))
    colnames(tmp) = c('c1', 'c2')
    expression_data = cbind(expression_data, tmp)
    
    comparisons = inner_join(comparisons, expression_data, by = c("c1", "c2"))
    test = function(){
        db_regions = as.character(comparisons[2,1])
        deseq_results = as.character(comparisons[2,4])
        combine_expression(db_regions, deseq_results)
    }


}

# }}}

# {{{
#setwd("/nfs/research2/bertone/user/mxenoph/hendrich/enhancers/hendrichChIP/macs/")

#SampleSheet must have SampleID,Condition,Treatment,Replicate,bamReads,bamControl,Peaks (last 3 are the path to the respective bed files)
DB<-function(sheet, fdr=0.05){
   library(DiffBind)

   output <- gsub(".csv", 'pdf', gsub(".*/", '', sheet))

   #Cause it will otherwise complain and terminate when you run it on a cluster
   pdf(output)
   #working object is of dba class

   #Check how many Conditions are there
   line <- unlist(strsplit(tolower(readLines(sheet[1], n=1)), ','))
   if(!"condition" %in% line) stop("Sample Sheet does not contain a Condition column")

   protein <- dba(sampleSheet= sheet)
   #Will print the correlation heatmap using occupancy(peak caller score) data
   plot(protein)

   sheet_table <-read.table(sheet, sep=',')
   conditions <- sort(levels(sheet_table$condition), decreasing=T)

   nCond <- match("condition", line)
   tmp <- rep('NULL', length(line))
   tmp[nCond] <- NA
   nCond <- sort(levels(read.table(sheet, colClasses=tmp, sep=",", head=T)[[1]]), decreasing=T)

   #calculates a binding matrix, scores based on read counts (affinity scores) rather than the previous one that plotted confidence scores.
   #DBA_SCORE_READS_MINUS => read count for interval-read count for interval in control..might be a better way of doing this
   #bScaleControl=TRUE => scale the control reads based on relative library size
   wobj <- dba.count(wobj, score=DBA_SCORE_READS_MINUS, bScaleControl=TRUE, bParallel=TRUE)
   btag <- TRUE

   if(sum(wobj$masks[[nCond[1]]]) < 2 | sum(wobj$masks[[nCond[2]]]) < 2){
      wobj <- dba.contrast(wobj, wobj$masks[[nCond[1]]], wobj$masks[[nCond[2]]], nCond[1], nCond[2])
      btag <- FALSE
   }
   else{
      #step is actually optional cause it's set up
      wobj <- dba.contrast(wobj, categories=DBA_CONDITION)
   }


   wobj <- dba.analyze(
                       wobj,
                       method= DBA_DESEQ,
                       bSubControl= TRUE,
                       #bFullLibrarySize= TRUE,
                       #SOS if no replicates this option should be set to FALSE
                       bTagwise= btag,
                       bParallel= TRUE)

   wobj.DB <- dba.report(wobj,
                         method=DBA_DESEQ,
                         #only sites with an absolute Fold value greater than equal to this will be included in the report
                         #fold= 2,
                         th= fdr,
                         initString= fdr,
                         file= paste(gsub(".csv", '', gsub(".*/", '', sheet)), '.diffbind', sep=''))
   wobj.df <- dba.report(wobj, method=DBA_DESEQ, th=fdr, DataType = DBA_DATA_FRAME)


   write.table(wobj.df, file=gsub(".csv", 'pdf', gsub(".*/", '', sheet)), sep="\t", quote=FALSE, row.names=FALSE)
   #plot(wobj)
   #dba.plotMA(wobj)
   #PCA based on the afinity data for all sites
   #dba.plotPCA(wobj, DBA_CONDITION)
   #PCA based on the affiity data for the DB sites
   #dba.plotPCA(wobj, contrast=1, th=.05)
   dev.off()
   return(list("db"=wobj.DB))
}

#myPeakFiles <- list.files(path=getwd(), full.name=T )
#samplesheets must point to modified macs peaks.bed with 4 columns instead of five
mySampleSheets <- list.files("..", full.name=T, pattern="csv")

#system.time(DB(mySampleSheets[1]), .01)


#Comment in and out depending what you are running
#Mi2b_2i <- DB(mySampleSheets[1], 0.05)
#Mi2b_Epi <- DB(mySampleSheets[2], 0.05)
#Mi2b_serum <- DB(mySampleSheets[3], 0.05)



#source("filewith get.annotation")
#mm9.annot<-get.annotation('may2012.archive.ensembl.org', 'mmusculus_gene_ensembl')




# }}}
