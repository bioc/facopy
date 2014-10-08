\name{plotPCA}
\alias{plotPCA}

\title{
facopy PCA
}
\description{
Plot a PCA of the samples based on their similarity given a certain set of alterations.
}
\usage{
plotPCA(fad, alteration, varName, sel = NULL, 
        design = c("binary", "versus", "lvog")[1], 
        do.plot = TRUE, by.size = TRUE, cex = 4)
}

\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data.
}
  \item{alteration}{
  A character describing the kinds of alteration to include. It should be one of the following:\cr
  - \code{amplifications }{All amplifications (CN>2).}\cr
  - \code{deletions }{All deletions (CN<2).}\cr
  - \code{loh }{All loss of heterozygosity (LOH), regardless of copy number.}\cr
  - \code{cnas }{All copy number alterations (CN<>2).}\cr
  - \code{any }{Any kind of alteration.}\cr
  - \code{all }{Any kind of alteration, same as \code{any}.}\cr
  - \code{onlygain }{Only non-LOH amplifications.}\cr
  - \code{someloss }{All deletions plus LOH alterations.}\cr
}
  \item{varName}{
A character indicating the variable of interest within the \code{facopyInfo} object. Points representing samples in the PCA will be colored according to the classification in such variable.\cr
Call \code{summary} on your \code{facopyInfo} object to see the names of defined variables.
}
  \item{sel}{
A character vector with the chromosome arms that will be taken into account in the analysis. Example: \code{c("1q","9p")}.
}
  \item{design}{
Depending on the chosen alteration, different designs are available. The simplest design is \code{binary}: an alteration exists or it does not. The \code{versus} design, for CNAs, assigns a value of -1, 0 or 1 depending on whether a deletion, no copy number change or an amplification exists for a given feature. The \code{vlog} design, for all (any) alterations, assigns a value of -1, 0 or 1 depending on whether a deletion or LOH, no copy number change or an amplification without LOH exists.
}
  \item{do.plot}{
A logical indicating whether to produce graphical output.
}
  \item{by.size}{
A logical indicating whether the width of the points' border represents the frequency of the selected alterations.
}
  \item{cex}{
A numerical value giving the amount by which points in the plot should be scaled relative to the base size. Default is 4. 
}
}
\value{
The results of the PCA are returned and the graphical output is optional.
}

\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

pca = plotPCA(myStudy, "any", "stage")
head(pca$eig)
}

\keyword{ plots }
