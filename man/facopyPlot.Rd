\name{facopyPlot}
\alias{facopyPlot}

\title{
Plot a Genome-Wide Overview of Association Results
}
\description{
The function plots alteration frequencies by genomic feature, broken down by variable value (if the variable is discrete). It also allows to plot relevant alteration frequencies within a specified external dataset.
}
\usage{
facopyPlot(fad, alteration, varName, db = NULL)
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
A character indicating the variable of interest within the \code{facopyInfo} object. Call \code{summary} on your \code{facopyInfo} object to see the names of defined variables.
}
  \item{db}{
An optional string representing the name of a database whose data will be overlaid in the output plot. Typically, the format is "[\code{database}]_[\code{dataset}]". The total amplification plus deletion frequencies will be displayed unless \code{alteration} indicates either amplification or deletion. In such cases, only the matching alterations are displayed. 
}
}
\details{
Alterations in the selected external database, if selected, are depicted as grey overlaid bars.
}
\value{
A plot is generated in the graphics device. Nothing is returned.
}
\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

facopyPlot(myStudy, "amp", "stage", db="gsk_colon")
}

\keyword{ plots }
