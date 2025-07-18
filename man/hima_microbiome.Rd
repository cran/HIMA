% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hima_microbiome.R
\name{hima_microbiome}
\alias{hima_microbiome}
\title{High-dimensional mediation analysis for compositional microbiome data}
\usage{
hima_microbiome(
  X,
  OTU,
  Y,
  COV = NULL,
  FDRcut = 0.05,
  verbose = FALSE,
  parallel = FALSE,
  ncore = 1
)
}
\arguments{
\item{X}{a vector of exposure. Do not use \code{data.frame} or \code{matrix}.}

\item{OTU}{a \code{data.frame} or \code{matrix} of high-dimensional Operational Taxonomic Unit (OTU) data (mediators).
Rows represent samples, columns represent variables.}

\item{Y}{a vector of continuous outcome. Binary outcome is not allowed. Do not use \code{data.frame} or \code{matrix}.}

\item{COV}{a \code{data.frame} or \code{matrix} of adjusting covariates. Rows represent samples, columns represent
microbiome variables. Can be \code{NULL}.}

\item{FDRcut}{Hommel FDR cutoff applied to select significant mediators. Default = \code{0.05}.}

\item{verbose}{logical. Should the function be verbose? Default = \code{FALSE}.}

\item{parallel}{logical. Enable parallel computing feature? Default = \code{FALSE}.}

\item{ncore}{number of cores to run parallel computing Valid when \code{parallel = TRUE}.}
}
\value{
A data.frame containing mediation testing results of significant mediators (FDR <\code{FDRcut}).
\describe{
    \item{Index: }{mediation name of selected significant mediator.}
    \item{alpha_hat: }{coefficient estimates of exposure (X) --> mediators (M) (adjusted for covariates).}
    \item{alpha_se: }{standard error for alpha.}
    \item{beta_hat: }{coefficient estimates of mediators (M) --> outcome (Y) (adjusted for covariates and exposure).}
    \item{beta_se: }{standard error for beta.}
    \item{IDE: }{mediation (indirect) effect, i.e., alpha*beta.}
    \item{rimp: }{relative importance of the mediator.}
    \item{pmax: }{joint raw p-value of selected significant mediator (based on Hommel FDR method).}
}
}
\description{
\code{hima_microbiome} is used to estimate and test high-dimensional mediation effects for compositional microbiome data.
}
\examples{
\dontrun{
# Note: In the following example, M1, M2, and M3 are true mediators.

data(MicrobiomeData)
pheno_data <- MicrobiomeData$PhenoData
mediator_data <- MicrobiomeData$Mediator

hima_microbiome.fit <- hima_microbiome(
  X = pheno_data$Treatment,
  Y = pheno_data$Outcome,
  OTU = mediator_data,
  COV = pheno_data[, c("Sex", "Age")],
  FDRcut = 0.05,
  verbose = TRUE
)
hima_microbiome.fit
}

}
\references{
1. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation effect selection in high-dimensional and compositional microbiome data.
Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955

2. Zhang H, Chen J, Li Z, Liu L. Testing for mediation effect with application to human microbiome data.
Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450
}
