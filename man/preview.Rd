\name{preview}
\alias{preview}
\title{
Alteration and Varible Preview
}
\description{
A wrapper function for \code{\link{variableSummary}}, \code{\link{alterationSummary}} and \code{\link{variableCor}}.
}
\usage{
preview(fad, folder = NULL)
}

\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data.
}
  \item{folder}{
Optional. A character specifying the name of the file to which to output the table with variable correlations.
}
}
\details{
A list of \code{\link{data.frame}} is always returned, and optionally written to files in a specified folder.
}
\value{
A list of \code{\link{data.frame}} with the results of calling the functions \code{\link{variableSummary}}, \code{\link{alterationSummary}} and \code{\link{variableCor}}.
}

\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

myCallsPreview = preview(myStudy)

myCallsPreview$byVar
head(myCallsPreview$byAlt)
myCallsPreview$varCor
}

\keyword{ overview }