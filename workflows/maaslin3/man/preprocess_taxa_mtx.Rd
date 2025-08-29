\name{preprocess_taxa_mtx}
\alias{preprocess_taxa_mtx}
\title{
Pre-process the taxa covariates for metatranscriptomics
}
\description{
Pre-process the taxa covariates for metatranscriptomics by
total-sum-scaling the taxa, matching the taxa to the RNAs coming from those
taxa, and then, for each sample in each feature:
\enumerate{
    \item{Log 2 transforming the taxon abundance if the taxon abundance is >=0}
    \item{Setting the taxon abundance to \code{log2([minimum non-zero relative
    abundance in the dataset] / 2)} if any of the corresponding RNA abundances
    are non-zero but the taxon abundance is zero}
    \item{Setting the taxon abundance to NA if both are zero}
}
When the taxon is present, the RNA data can be modeled as usual in MaAsLin
3 with
\code{log2(taxon)} as a covariate. When the taxon is not present, if any of
its RNA is present, we assume the taxon was missed due to finite read depth, 
so the taxon abundance is imputed with a small pseudo-count. When neither the
taxon nor RNA is
present, we assume the gene/microbe was not in the sample and therefore no
information about the transcription level can be obtained. Setting the taxon
covariate to NA has the effect of dropping the sample when fitting the relevant
feature's model in MaAsLin 3. Unlike most MaAsLin functions that will infer the
samples from the row names and column names, \bold{the \code{rna_table} must be
formated as samples (rows) by features (columns)}.
}
\usage{
preprocess_taxa_mtx(taxa_table, rna_table, rna_per_taxon)
}
\arguments{
\item{taxa_table}{The samples (rows) by features (columns) data frame of
taxon abundances to preprocess. These can be relative abundances or
counts.}
\item{rna_table}{The samples (rows) by features (columns) data frame of
RNA to preprocess. These can be relative abundances or counts.}
\item{rna_per_taxon}{A dataframe with the columns `RNA` and `taxon` with
one row per `RNA` column found in `rna_table` giving both the `RNA` column 
and which `taxon` column it corresponds to in `taxa_table`.}
}
\value{
A list containing the following named items:
    \enumerate{
    \item{\code{dna_table}: The table of log2 transformed taxon relative
    abundances with NAs for any feature-sample pairs for which both the taxon
    and RNA abundances were 0.}
    \item{\code{rna_table}: The table of total sum scaled RNA abundances.
    These are not log2 transformed.}
    }
}
\author{
    William Nickols<willnickols@g.harvard.edu>,\cr
    Jacob Nearing<nearing@broadinstitute.org>,\cr
    Maintainers: Lauren McIver<lauren.j.mciver@gmail.com>,\cr
}
\examples{
    taxa_in <- data.frame('tax1' = c(1, 2, 0, 4, 5), 
                    'tax2' = c(2, 3, 4, 5, 6), check.names = FALSE)
    rownames(taxa_in) <- paste0("sample", c(1:5))
    
    mtx_in <- data.frame('a' = c(1, 2, 3, 4, 5), 
                        'b' = c(2, 3, 4, 5, 0),
                        'c' = c(3, 4, 5, 6, 0), check.names = FALSE)
    rownames(mtx_in) <- paste0("sample", c(1:5))
    
    rna_per_taxon <- data.frame(RNA = c('a', 'b', 'c'),
                                taxon = c('tax1', 'tax1', 'tax2'))
    
    preprocess_out <- preprocess_taxa_mtx(taxa_in, mtx_in, rna_per_taxon)
}
