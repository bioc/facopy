\name{addVariables}
\alias{addVariables}

\title{
Add Variable Annotation
}
\description{
Attaches phenotypic information to a \code{facopyInfo} object that had been created from copy number data.
}
\usage{
addVariables(fad, varInfo, varTypes, 
             varColumns = NULL, varValues = NULL,
             varColumnsNames = NULL, varValuesNames = NULL,
             ...)
}
\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data. 
}
  \item{varInfo}{
A character or a \code{\link{data.frame}}. If it is a character, it should indicate the name of the file with information on the variables. If it is a code{\link{data.frame}}, a column named \code{code} should contain the samples' unique identifiers.
}
  \item{varTypes}{
A vector of characters with the same length as the number of variables whose information will be loaded. Variables can be either categorical or quantitative. Planning ahead for subsequent analyses, it might be best to regard ordered variables as categorical. Vector elements should equal either \code{"categorical"} or \code{"quantitative"}.
}
  \item{varColumns}{
An optional character vector. If the table with variable information contains extra columns, this parameter allows the specification of those columns with information on the variables to use.
}
  \item{varValues}{
An optional list of vectors. Each vector is either null or contains the subset of values to be used in subsequent analyses for the corresponding variable.
}
  \item{varColumnsNames}{
An optional character vector with the names of the variables, if different from the corresponding column names in the table. 
}
  \item{varValuesNames}{
An optional list of vectors with the names of the variable values, if different from the corresponding values in the columns. If a vector is null, the corresponding value names are not changed, which is useful for quantitative variables with many possible values.
}
  \item{\dots}{
Extra parameters passed to the function are in turn passed to \code{\link{read.table}}, such as the presence of header and the column separation character. Only used if \code{varInfo} is a character that indicates the name of the file with variable information.
}
}
\value{
\code{facopyInfo} object with incorporated phenotypic information. 
}
\author{
David Mosen-Ansorena
}
\examples{
data(myCalls)
data(myVariables)

myStudy = addVariables(myCalls, myVariables, c("continuous","categorical")) 
}

\keyword{ preprocessing }
