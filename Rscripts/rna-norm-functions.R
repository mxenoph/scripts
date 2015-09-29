# both counts and effective length are dataframes of ncol = samples and 1 respectively
counts_to_TPM = function(counts, effective_length){
    rate = log(counts) - log(effective_length)
    denom = sapply(1:ncol(rate), function(x) log((sum(exp(rate[[x]])))))
#    denom = log((sum(exp(rate))))
    exp(rate - denom + log(1e6))
}

counts_to_FPKM = function(counts, effective_length){
    N = sum(counts)
    exp( log(counts) + log(1e9) - log(effective_length) - log(N) )
}

FPKM_to_TPM = function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

counts_to_EffCounts = function(counts, len, effective_length){
    counts * (len / effective_length)
}

compute_effective_length = function(gene_length, fragment_length){
    x = gene_length %>% mutate(efflength = length - fragment_length +1)
    return(x)
}

get_effective_length = function(gene, fragment_length){
    effective_length = gene - fragment_length + 1
}
