\name{alterationSummary}
\alias{alterationSummary}

\title{
Alteration Summary
}
\description{
Summary of copy number calls in a \code{facopyInfo} object, by alteration type and chromosome arm.
}
\usage{
alterationSummary(fad, filename = NULL)
}
\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data.
}
  \item{filename}{
Optional. A character specifying the name of the file to which to output the alteration summary table.
}
}
\details{
A \code{\link{data.frame}} is always returned, and optionally written to a file.
}
\value{
A \code{\link{data.frame}} with alteration frequencies by arm and alteration type.
}
\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

alterationSummary(myStudy)
}

\keyword{ overview }
