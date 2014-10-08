\name{facopy}
\alias{facopy}

\title{
Feature-Based Copy Number Association Analysis
}
\description{
Main function in the facopy package. It performs statistical association between copy number data and further variables at each genomic feature of interest.
}
\usage{
facopy(fad, alteration, model, nullModel = NULL,
       modelPart = c("response","predictor","unknown","whole")[1],
       strata = NULL, toOrdered = NULL, toIntervals = NULL, 
       sel = NULL, plot = FALSE, pvalThr = 0.05, db = NULL, 
       link = c("logit", "probit")[1], parametric = FALSE,
       design = c("binary", "versus", "lvog")[1],
       FUN, ...)
}

\arguments{
  \item{fad}{
\code{facopyInfo} object with a certain study's \code{facopy} data. Feature information should have been added beforehand on this object (see \code{\link{addFeatures}}).
}
  \item{alteration}{
A character, the name of the combination of alterations to be considered for the analysis. \cr
It should be one of the following:\cr
  - \code{amplifications }{All amplifications (CN>2).}\cr
  - \code{deletions }{All deletions (CN<2).}\cr
  - \code{loh }{All loss of heterozygosity (LOH), regardless of copy number.}\cr
  - \code{cnas }{All copy number alterations (CN<>2).}\cr
  - \code{any }{Any kind of alteration.}\cr
  - \code{all }{Any kind of alteration, same as \code{any}.}\cr
  - \code{onlygain }{Only non-LOH amplifications.}\cr
  - \code{someloss }{All deletions plus LOH alterations.}\cr
}
  \item{model}{
A character. Model, or part of it (response or predictor), whose association will be measured.\cr
If \code{modelPart="response"}, the name of a single variable, representing the response of the association model.\cr
If \code{modelPart="predictor"}, a linear predictor of the copy number, using a combination of the variables.\cr
If \code{modelPart="unknown"}, the association is measured as the strenght of the relationship between copy number and the specified combination of variables. The model is limited to \code{x1 + ... + xn | strata}, see \code{\link{independence_test}}.\cr
If \code{modelPart="whole"}, a character representation of a formula in compact symbolic form. Use the @ (at) symbol to refer to the copy number variable.\cr
}
  \item{nullModel}{
A character. The null model against which to evaluate the fitness of the association model. Only used if \code{modelPart="predictor"} or \code{modelPart="whole"}.
}
  \item{modelPart}{
A character. Indicates what part of the association model is defined in the \code{model} parameter.
}
  \item{strata}{
A character. The name of a categorical variable that defines stratification blocks in the association model. Only used if \code{modelPart="response"} or \code{modelPart="unknown"}.
}
  \item{toOrdered}{
Certain categorical variables can also be understood as ordered. This parameter takes a list of named vectors, where the name of each vector is a variable name and its contents reflect the quantification of the variable values, in the same order as defined in \code{\link{addVariables}}.
}
  \item{toIntervals}{
Quantitative variables can be broken down into intervals. This parameter takes a list of named vectors, where the name of each vector is a variable name and its contents reflect the breaks of the variable values (excluding bottom and top limits).
}
  \item{sel}{
A vector of feature names on which to perform the association. Leave to \code{NULL} for genome-wide association over all the features specified in \code{\link{addFeatures}}. Plotting will only be done if \code{sel=NULL}.
}
  \item{plot}{
A logicil indicating whether to output a composite plot with an arm-wise display of genome-wide alteration frequencies. If the model consists of a single variable, frequencies will be broken down by variable value. Features with significant associations and additional information pulled from external databases can be displayed as overlaid layers.
}
  \item{pvalThr}{
Significant associations under this threshold will be shaded in the output plot. Only used if \code{plot=TRUE}.
}
  \item{db}{
An optional string representing the name of a database whose data will be overlaid in the output plot. Typically, the format is "[\code{database}]_[\code{dataset}]". The total amplification plus deletion frequencies will be displayed unless \code{alteration} indicates either amplification or deletion. In such cases, only the matching alterations are displayed. Only used if \code{plot=TRUE}. See \code{\link{getFacopyInfo}} for a list of available data sets.
}
  \item{link}{
A character, link function to be used with the multinomial error distribution in logistic regression models. See \code{\link{glm}}.
}
  \item{parametric}{
A logical that indicates whether to perform one-way ANOVA instead of Kruskal-Wallis in the association of copy number with quantitative variables.
}
  \item{design}{
Depending on the chosen alteration, different designs are available. The simplest design is \code{binary}: an alteration exists or it does not. The \code{versus} design, for CNAs, assigns a value of -1, 0 or 1 depending on whether a deletion, no copy number change or an amplification exists for a given feature. The \code{vlog} design, for all (any) alterations, assigns a value of -1, 0 or 1 depending on whether a deletion or LOH, no copy number change or an amplification without LOH exists.  
}
  \item{FUN}{
A function that tests a model. Only used if \code{modelPart="whole"}. Functions from the \code{coin} package and those that inherit from either the \code{lm} (such as \code{glm}) are supported, as well as those that directly return a pvalue. Thus, the function of interest can be wrapped in a wrapper function that provides the pvalue.
}
  \item{\dots}{
Further arguments for the \code{FUN} function. Only used if \code{modelPart="whole"}.
}
}
\details{
Only the \code{facopyInfo} object, the alteration type and a simple model (e.g. the name of a variable) are required. The rest of the parameters tune up the association model and control the graphical output.\cr
Alterations in the selected external database, if selected, are depicted as grey overlaid bars. Significant regions are depicted in turn as overlaid rectangles that go from top to bottom.
}
\value{
A \code{\link{data.frame}} with the following columns:
\item{feature }{Name of the genomic feature.}
\item{p_value }{Pvalue from the association test under the given model at the genomic feature.}
\item{chr_q_arm }{Chromosome and arm in which the genomic feature lies.}
\item{bp_st }{Starting genomic position of the feature within the chromosome.}
\item{bp_en }{Ending genomic position of the feature within the chromosome.}
}

\author{
David Mosen-Ansorena
}

\examples{
data(myStudy) # load example study

genes = facopy(myStudy, "amp", "stage")
head(genes)
}
\keyword{ association }
