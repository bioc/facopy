\name{variableCor}
\alias{variableCor}

\title{
Variable Correlations
}
\description{
Correlations between pairs of variables in a \code{facopyInfo} object.
}
\usage{
variableCor(fad, filename = NULL)
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
A \code{\link{data.frame}} with alteration variable correlations, with a combination of the appropriate parametric and non-parametric tests.
}
\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

variableCor(myStudy)
}
\keyword{ overview }
