\name{addFeatures}
\alias{addFeatures}

\title{
Add Feature Annotation
}
\description{
Adds feature annotation to a \code{facopyInfo} object and calculates the overlapping frequency in the sample set for every combination of copy number and feature.
}
\usage{
addFeatures(fad, what = c("ensembl", "cancergene", "oncogene", 
                          "tumorsupressor", "lincRNA", "mirnas")[1], 
            genome = c("hg18", "hg19", "mm8")[1], 
            lMargin = 0, rMargin = 0, minoverlap = 1, 
            data = NULL)
}

\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data. Variable information should have been added beforehand on this object (see \code{\link{addVariables}}).
}
  \item{what}{
A character indicating the name of the feature set to use in subsequent analyses. If the name belongs to a set of features integrated in \code{facopy}, annotation is automatically loaded and the parameter \code{data} is not used.
}
  \item{genome}{
Reference genome build used in the copy number calling step prior to running \code{facopy}. Some feature sets may not be available depending on the selected genome.
}
  \item{lMargin}{
Number of base pairs with which to expand the left flanking side of each feature, increasing its chance to overlap alterations.
}
  \item{rMargin}{
Number of base pairs with which to expand the right flanking side of each feature, increasing its chance to overlap alterations.
}
  \item{minoverlap}{
Minimum overlap, in base pairs, between each feature and alteration.
}
  \item{data}{
If the genomic features of interest are not integrated within \code{facopy}, their information can be manually loaded using this parameter.
Either a character or a \code{\link{data.frame}}. If a character, it indicates the file with the information, with headers and tab-delimited. In either case, the information shoul follow a structure with the following columns:\cr
- \code{chr }{Chromosome in which the genomic feature lies.}\cr
- \code{bp_st }{Starting genomic position of the feature within the arm.}\cr
- \code{bp_en }{Ending genomic position of the feature within the arm.}\cr
- \code{feature }{Name of the feature.}\cr
- \code{chr_q_arm }{Chromosome arm in which the feature lies.}
}
}
\value{
\code{facopyInfo} object with information on both variable and feature annotations. 
}
\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

myStudy = addFeatures(myStudy, "oncogene", "hg18")
}

\keyword{ preprocessing }
