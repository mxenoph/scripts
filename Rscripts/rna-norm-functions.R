library = function (...) suppressMessages(base::library(...))

#' Calculates the total length of all blocks in features (i.e. genes or transcripts)# {{{
#' @param gtf path to an ensembl gtf file
#' @param feature the feature type for which effective lenght is calculated (gene or transcript)
#' works best for defaults the gtf produced from get-annotation-rscript.r
get_lengths = function(gtf, feature = 'ensembl', type = 'gff3') {
    x = c('GenomicRanges', 'rtracklayer', 'tools', 'BiocParallel')
    lapply(x, library, character.only=T)
    output = file.path(dirname(gtf),
                       paste0(basename(file_path_sans_ext(gtf)),
                              '.', feature, '-lengths.tsv'))

    gtf = import(gtf, format = type)
    per_feature = split(gtf, as.factor(values(gtf)[[paste0(feature, '_id')]]))

    # R is running interactively so parallel::detectCores() will count the cores on# {{{
    # the machine the job is running. In that case just use 4 cores
    if (interactive()) {
        ncores = 4
    } else {
        ncores = parallel::detectCores()
    }
    register(MulticoreParam(workers = ncores))
    # }}}

    #per_feature = bplapply(per_feature, reduce, BPPARAM = bpparam())
    feature_length = sapply(per_feature, function(x) sum(width(x)))

    write.table(as.data.frame(feature_length),
                file = output, col.names = FALSE, sep = "\t", quote = FALSE)
    return(feature_length)
}
# }}}

get.lengths = get_lengths

#' Compute TPMs from normalised ounts# {{{
#' @param counts can be a vector/matrix/dataframe
#' @param effective_length is a vector of length length(counts) or nrow(counts)
counts_to_TPM = function(counts, effective_length){
    library(dplyr)
    # not a matrix provided
    if (is.null(dim(counts))){
        if(all.equal(names(counts), names(effective_length)) != TRUE)
           stop('Counts and effective lengths are not calculated on the same features')

        rate = log(counts / effective_length)
        denom = log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
    } else {
        if(all.equal(rownames(counts), names(effective_length)) != TRUE)
           stop('Counts and effective lengths are not calculated on the same features')

        tpms = as.data.frame(counts) %>% 
                mutate_each(funs(exp((log(. / effective_length)) - log(sum(exp(log(. / effective_length)))) + log(1e6))))
        rownames(tpms) = rownames(counts)
        tpms
    }
}# }}}

# TODO: find scripts using this function and use counts_to_FPKM instead
# as it has more checks and uses dplyr instead of for loops
compute.FPKMS=function(counts, length){# {{{
    #Give it a matrix for exons or introns and a matrix with the combined lengths of the respective features.
    total = NULL
    for(c in 1:ncol(counts)){
        total = cbind(total, sum(counts[,c]))
    }
    total = total/10^6

    fpkm = merge(counts, length, by= "row.names")
    rownames(fpkm) = fpkm[,1]
    remove = "Row.names"
    fpkm = fpkm[, !colnames(fpkm) %in% remove]
    myncol = ncol(fpkm)
    m = (ncol(fpkm))-1
    for(i in 1:m){
        fpkm[,i] = fpkm[,i] / ((fpkm[,myncol] / 10^3) * (total[1,i]))
    }
    return(list(matrix = fpkm, total = total))
}# }}}

#z' Compute FPKMs from raw counts# {{{
#' @param counts can be a vector/matrix/dataframe
#' @param effective_length is a vector of length length(counts) or nrow(counts)
counts_to_FPKM = function(counts, effective_length){
    library(dplyr)

    if (is.null(dim(counts))){
        # If not the same order or not of the same length this will return a string 
        # with the number of mismatches
        if(all.equal(names(counts), names(effective_length)) != TRUE)
           stop('Counts and effective lengths are not calculated on the same features')

        N = sum(counts)
        exp( log(counts) + log(1e9) - log(effective_length) - log(N) )
    } else {
        if(all.equal(rownames(counts), names(effective_length)) != TRUE)
           stop('Counts and effective lengths are not calculated on the same features')

        fpkms = as.data.frame(counts) %>%
                mutate_each(funs(exp( log(.) + log(1e9) - log(effective_length) - log(sum(.)) )))
        rownames(fpkms) = rownames(counts)
        fpkms
    }
}#}}}

FPKM_to_TPM = function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

counts_to_EffCounts = function(counts, len, effective_length){
    counts * (len / effective_length)
}

compute_effective_length = function(feature_length, fragment_length){
    library(dplyr)
    # Useful discussion in https://groups.google.com/forum/embed/#!topic/rsem-users/gDh0OBLK6_4
    # After discussing with Mitra the effective length for *the gene* should be coomputed as the 
    # mean(transcript size) - mean(fragment size). This function expects that lengths is the mean(transcript size)
    x = feature_length %>% mutate(efflength = lengths - fragment_length +1)
    # Dealing with genes that have sum(exons) < fragment size
    # These can only give on fragment in library prep
    x = x %>% mutate(efflength = ifelse(efflength <= 0, 1, efflength))
    rownames(x) = rownames(feature_length)
    return(x)
}

