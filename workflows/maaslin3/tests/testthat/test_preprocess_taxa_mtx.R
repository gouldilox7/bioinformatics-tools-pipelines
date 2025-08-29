library(testthat)
library(maaslin3)

taxa_in <- data.frame('tax 1' = c(1, 0, 0, 4, 5), 
                    'tax 2' = c(2, 3, 4, 5, 6), check.names = FALSE)
rownames(taxa_in) <- paste0("sample", c(1:5))

mtx_in <- data.frame('a' = c(1, 0, 0, 4, 5), 
                    'b' = c(2, 3, 0, 5, 0),
                    'c test' = c(3, 4, 5, 6, 0), check.names = FALSE)
rownames(mtx_in) <- paste0("sample", c(1:5))

rna_per_taxon <- data.frame(RNA = c('a', 'b', 'c test'),
                            taxon = c('tax 1', 'tax 1', 'tax 2'))

data_in_tss <- data.frame(t(apply(taxa_in, MARGIN = 1, 
                                FUN = function(x){x / sum(x)})), 
                        check.names = F)
data_in_tss <- data_in_tss[,c('tax 1', 'tax 1', 'tax 2')]
colnames(data_in_tss) <- c('a', 'b', 'c test')
data_in_tss[3, 1:2] <- NA # Because all RNA for taxa 1 are missing in sample 3
data_in_tss[data_in_tss == 0] <- min(data_in_tss[data_in_tss > 0], na.rm=T) / 2
data_in_tss_log <- log2(data_in_tss)

expect_that(preprocess_taxa_mtx(taxa_in, mtx_in, rna_per_taxon)$dna_table, 
            equals(data_in_tss_log))

data_in_tss <- data.frame(t(apply(mtx_in, MARGIN = 1, 
                            FUN = function(x){x / sum(x)})), check.names = FALSE)

expect_that(preprocess_taxa_mtx(taxa_in, mtx_in, rna_per_taxon)$rna_table, 
                                equals(data_in_tss))



