\name{plotBar}
\alias{plotBar}

\title{
Overview of Alterations by Chromosome Arm
}
\description{
The function plots amplification and deletion frequencies by chromosome arm, as well as LOH frequencies
}
\usage{
plotBar(fad, byFeature = TRUE, sel = NULL, selColors = NULL, 
        ylim = c(-1, 1), baseColor = "black", 
        varName = NULL, value = NULL)
}

\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data.
}
  \item{byFeature}{
A logical. If \code{TRUE}, arm frequencies are calculated as the mean alteration frequencies across genomic features in each arm. If \code{TRUE}, they are calculated as the mean proportion of altered base pairs in each arm.
}
  \item{sel}{
A character vector with the chromosome arms that will be highlighted using the colors in \code{selColors}. Example: \code{c("1q","9p")}.
}
  \item{selColors}{
A vector of colors with the same length as the parameter \code{sel}. Each color in the vector will be used for the selected chromosome arm in the same position within \code{sel}.
}
  \item{ylim}{
A two integer vector indicating frequency limits in the plot. Deletion frequency is indicated with a negative number. Default is \code{c(-1,1)}.
}
  \item{baseColor}{
Base color of the triangles that indicate frequencies in the plot. Default is \code{"black"}.
}
  \item{varName}{
Either \code{NULL} or a character indicating the variable of interest within the \code{facopyInfo} object. Call \code{summary} on your \code{facopyInfo} object to see the names of defined variables.\cr
If \code{NULL}, the frequencies across all samples are plotted. Otherwise, only samples with the variable value in the parameter \code{value} are selected for the computation. In this latter case, only categorical variables are allowed.
}
  \item{value}{
Value that the samples should have for the selected variable. Only these samples are selected for the frequency calculation.
}
}
\details{
Amplification frequencies are depicted by chromosome arm over the horizontal line, while deletion frequencies appear below it. The frequencies of alterations with LOH are depicted as smaller triangles.
}
\author{
David Mosen-Ansorena
}
\examples{
data(myStudy) # load example study

# select some chromosome arms to highlight
myArms = c("8q","13q","20q","8p","18q")
myColors = c(rainbow(15)[1:3], rainbow(15)[10:11])

plotBar(myStudy, TRUE, myArms, myColors, ylim=c(-0.5,1))
}

\keyword{ plots }
