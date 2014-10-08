\name{variableSummary}
\alias{variableSummary}
\title{
Variable Summary
}
\description{
For each pair of values in a variable, it performs a statistical test to see whether samples with such variable values present significantly different spans of a certain alteration type. 
}
\usage{
variableSummary(fad, filename = NULL)
}

\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data.
}
  \item{filename}{
Optional. A character specifying the name of the file to which to output the table with variable correlations.
}
}
\details{
A \code{\link{data.frame}} is always returned, and optionally written to a file.
}
\value{
A \code{\link{data.frame}} with alteration spans and significance for each pair of values in a variable and alteration type.
}

\author{
David Mosen-Ansorena
}
\examples{
data(myStudy) # load example study

variableSummary(myStudy)
}

\keyword{ overview }