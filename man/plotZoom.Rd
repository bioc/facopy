\name{plotZoom}
\alias{plotZoom}
\title{
Zoom In to Plot a Chromosome Arm or a Genomic Feature
}
\description{
The function shows the selected alterations either in a chromosome arm or near a genomic feature, color-coded by variable value.
}
\usage{
plotZoom(fad, what = c("feature", "arm"), name, 
         alteration, varName, margin = 0)
}

\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data.
}
  \item{what}{
Charcter indicating whether to plot selected alterations either near a genomic feature or in a chromosome arm. Default is \code{feature}.
}
  \item{name}{
Name of the genomic feature or chromosome arm of interest. Examples: \code{"MIRN181A2"}, \code{"14q"}.
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
A character indicating the variable of interest within the \code{facopyInfo} object. Bars representing alterations will be colored according to the classification of the corresponding samples in such variable.\cr
Call \code{summary} on your \code{facopyInfo} object to see the names of defined variables.
}
  \item{margin}{
If parameter \code{what} is \code{"feature"}, a margin downstream and upstream of the genomic feature to consider when searching for overlapping alterations.
}
}
\details{
One line per sample in the study is displayed in the case of chromosome arms, whereas only those samples with alterations are displayed in the case of genomic features.
}

\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

plotZoom(myStudy, "feat", "RAD51L1", "amp", "stage")

plotZoom(myStudy, "arm", "8p", "del", "stage")
}

\keyword{ plots }
