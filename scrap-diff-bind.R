
# if working with meryem files just give one of the files (not the -de.tsv though)# {{{
retrieve_annotation = function(fs){

    # If de files from Meryem, gene coordinates are in deseq files. # {{{
    # No need to load get_annotation(assembly) to get gene coordinates
    # chromosome names and strands needs to be changed
    if(!missing(fs)){
        fs = read.delim(fs, head = T, sep = "\t", stringsAsFactors = F)
        source("~/source/Rscripts/granges-functions.R")
        genes = with((fs %>% dplyr::select(chromosome_name:end_position)),
                     GRanges(seqnames = chromosome_name,
                             IRanges(start_position, end_position),
                             strand = strand))
        values(genes) = (fs %>% dplyr::select(starts_with("gene")))
        names(genes) = genes$gene_id
        genes = change_seqnames(genes)
    } else {
        #ToDO: load genes from get_annotation(assembly)
    }
    # }}}
    return(genes)
}
# }}}
aggregate_count_data = function(fs, type = "binding"){ # {{{
    flag = FALSE
    for (x in 1:length(fs)) {
        tmp = read.delim(fs[x], head = T, sep = "\t", stringsAsFactors = F)

        rename_map = c(paste("log2FoldChange", names(fs)[x], sep = "."),
                       paste("FDR", names(fs)[x], sep = "."))
        combine_by = c('seqnames', 'strand', 'start', 'end')

        if(type == "expression") {# {{{
            # Don't try to read it in and select on one go as it won't work
            tmp = tmp %>% dplyr::select(gene_id, chromosome_name:end_position,
                                        log2FoldChange, padj:ncol(.)) %>%
                  rename_(.dots = setNames(list("log2FoldChange", "padj",
                                                "chromosome_name", "start_position", "end_position"),
                          c(rename_map, "seqnames", "start", "end")))

            combine_by = c('gene_id', combine_by)
        } else {
            tmp = tmp %>% dplyr::select(seqnames:end, strand, Fold, FDR:ncol(.)) %>%
                  rename_(.dots = setNames(list("Fold", "FDR"), rename_map))
        }
        # }}}

        if(flag) {# {{{
            ofs = inner_join(ofs, tmp, by = combine_by)
            # since samples at h0 will be in all comparisons expression columns for those will be
            # duplicated and dplyr add '.x' or '.y' to them depending from which df they came for
            # remove extension, find duplicates and remove. This would work nicer with rename() but I
            # can't get it to work on selected columns
            colnames(ofs) = gsub("\\.x|\\.y$", "", colnames(ofs))
            discard = ofs %>% names(.) %>% duplicated()
            ofs = ofs[,!discard]
        } else {
            ofs = tmp
            flag = TRUE
        }
    }# }}}
    return(ofs)
}
# }}}
