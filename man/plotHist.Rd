\name{plotHist}
\alias{plotHist}
\title{
Stacked Histograms
}
\description{
The function plots a stacked histogram of chromosome arm alterations by variable value
}
\usage{
plotHist(fad, alteration, varName, 
         sel = NULL, selColors = NULL, selOnly = FALSE, 
         baseColor = "black", bin = 0.05, xmax = 1, ymax)
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
A character indicating the variable of interest within the \code{facopyInfo} object. A stacked histogram will be generated for every value of the variable. Only available for discrete variables.\cr
Call \code{summary} on your \code{facopyInfo} object to see the names of defined variables.
}
  \item{sel}{
A character vector with the chromosome arms that will be highlighted using the colors in \code{selColors}. Example: \code{c("1q","9p")}.
}
  \item{selColors}{
A vector of colors with the same length as the parameter \code{sel}. Each color in the vector will be used for the selected chromosome arm in the same position within \code{sel}.
}
  \item{selOnly}{
A logical indicating whether to only plot the frequencies of selected chromosome arms.
}
  \item{baseColor}{
Base color of the triangles that indicate frequencies in the plot. Default is \code{"black"}.
}
  \item{bin}{
The width of the histogram cells. Default is 0.05.
}
  \item{xmax}{
The maximum limit of the x axis in the plot. Default is 1.
}
  \item{ymax}{
The maximum limit of the y axis in the plot.
}
}
\details{
Alteration frequencies are calculated genomic feature-wise within each chromosome arm.
}
\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

# select some chromosome arms to highlight
myArms = c("8q","13q","20q","8p","18q")
myColors = c(rainbow(15)[1:3], rainbow(15)[10:11])

plotHist(myStudy, "amp", "stage", myArms, myColors, bin=0.1, ymax=80)
}


\keyword{ plots }
