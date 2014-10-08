\name{calculateCor}
\alias{calculateCor}

\title{
Calculate Correlation Between Copy Number and Expression
}
\description{
Calculates the gene-wise correlation between copy number and expression data, which might come from: (1) the \code{\link{facopyInfo}} object and user-provided expression data or (2) from an external database.
}
\usage{
calculateCor(fad, exprProfile, db = NULL)
}
\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data. 
}
  \item{exprProfile}{
Either a character or a \code{\link{data.frame}}. If \code{db} is \code{NULL}, a \code{\link{data.frame}} with the expression data or a character that indicates the file with the expression data, with headers and tab-delimited. Otherwise, a character indicating the the name of the expression profile in the Cancer Genomics Data Server. See \code{http://www.cbioportal.org/public-portal/cgds_r.jsp}.
}
  \item{db}{
A character indicating the the name of the dataset, in the Cancer Genomics Data Server, from which to get the expression profile. See \code{http://www.cbioportal.org/public-portal/cgds_r.jsp}.
}
}
\details{
Only use this function and \code{\link{facopyEnrichment}} if you selected some kind of gene collection as genomic features.
}
\value{
The returned object is used to select genes in the enrichment process encapsulated in the \code{\link{facopyEnrichment}} function.
}
\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

eCor = calculateCor(myStudy, "mrna_merged_median_Zscores", "coadread_tcga_pub")
head(eCor$cor)
}

\keyword{ enrichment }
