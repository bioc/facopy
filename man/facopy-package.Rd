\name{facopy-package}
\alias{facopy-package}
\docType{package}
\title{
facopy: feature-based association and gene-set enrichment for copy number alteration analysis in cancer
}
\description{
facopy is an R package for fine-tuned cancer CNA association modeling. Association is measured directly at the genomic features of interest and, in the case of genes, downstream gene-set enrichment analysis can be performed thanks to novel internal processing of the data. The software opens a way to systematically scrutinize the differences in CNA distribution across tumoral phenotypes, such as those that relate to tumor type, location and progression. Currently, the output format from 11 different methods that analyze data from whole-genome/exome sequencing and SNP microarrays, is supported. Multiple genomes, alteration types and variable types are also supported.
}
\details{
\tabular{ll}{
Package: \tab facopy\cr
Type: \tab Package\cr
License: \tab CC BY-NC 4.0\cr
}
To analyze your data with facopy, first read copy number data through \code{\link{readCNData}}. Then, \code{\link{addVariables}} and \code{\link{addFeatures}}, and you are ready to start scrutinizing the data in multiple ways thanks to the many provided functions.
}
\author{
David Mosen-Ansorena
}

\keyword{ package }

