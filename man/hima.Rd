% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hima.R
\name{hima}
\alias{hima}
\title{High-dimensional Mediation Analysis}
\usage{
hima(
  formula,
  data.pheno,
  data.M,
  mediator.type = c("gaussian", "negbin", "compositional"),
  penalty = c("DBlasso", "MCP", "SCAD", "lasso"),
  quantile = FALSE,
  efficient = FALSE,
  scale = TRUE,
  sigcut = 0.05,
  contrast = NULL,
  subset = NULL,
  verbose = FALSE,
  parallel = FALSE,
  ncore = 1,
  ...
)
}
\arguments{
\item{formula}{an object of class \code{formula} representing the overall effect model to be fitted, specified as \code{outcome ~ exposure + covariates}.
The "exposure" variable (the variable of interest) must be listed first on the right-hand side of the formula.
For survival outcomes specified using \code{Surv()}, the exposure should be the first
variable after the \code{\link[base:tilde]{~}}.}

\item{data.pheno}{a data frame containing the exposure, outcome, and covariates specified in the formula. Variable names in \code{data.pheno} must match those 
in the formula. When \code{scale = TRUE}, the exposure and covariates will be scaled (the outcome retains its original scale).}

\item{data.M}{a \code{data.frame} or \code{matrix} of high-dimensional mediators, with rows representing samples and columns representing mediator variables. 
When \code{scale = TRUE}, \code{data.M} will be scaled.}

\item{mediator.type}{a character string indicating the data type of the high-dimensional mediators (\code{data.M}). Options are: \code{'gaussian'} (default): 
for continuous mediators. \code{'negbin'}: for count data mediators modeled using the negative binomial distribution (e.g., RNA-seq data). \code{'compositional'}: 
for compositional data mediators (e.g., microbiome data).}

\item{penalty}{a character string specifying the penalty method to apply in the model. Options are: \code{'DBlasso'}: De-biased LASSO (default). \code{'MCP'}: 
Minimax Concave Penalty. \code{'SCAD'}: Smoothly Clipped Absolute Deviation. \code{'lasso'}: Least Absolute Shrinkage and Selection Operator. Note: Survival HIMA and microbiome 
HIMA can only be performed with \code{'DBlasso'}. Quantile HIMA and efficient HIMA cannot use \code{'DBlasso'}; they always apply \code{'MCP'}.}

\item{quantile}{logical. Indicates whether to use quantile HIMA (\code{hima_quantile}). Default is \code{FALSE}. Applicable only for classic HIMA with a continuous outcome and 
\code{mediator.type = 'gaussian'}. If \code{TRUE}, specify the desired quantile(s) using the \code{tau} parameter; otherwise, the default \code{tau = 0.5} (i.e., median) is used.}

\item{efficient}{logical. Indicates whether to use efficient HIMA (\code{hima_efficient}). Default is \code{FALSE}. Applicable only for classic HIMA with a continuous outcome and 
\code{mediator.type = 'gaussian'}.}

\item{scale}{logical. Determines whether the function scales the data (exposure, mediators, and covariates). Default is \code{TRUE}. Note: For simulation studies, set 
\code{scale = FALSE} to avoid estimate compression (i.e., shrinkage of estimates toward zero due to scaling).}

\item{sigcut}{numeric. The significance cutoff for selecting mediators. Default is \code{0.05}.}

\item{contrast}{a named list of contrasts to be applied to factor variables in the covariates (cannot be the variable of interest).}

\item{subset}{an optional vector specifying a subset of observations to use in the analysis.}

\item{verbose}{logical. Determines whether the function displays progress messages. Default is \code{FALSE}.}

\item{parallel}{logical. Enable parallel computing feature? Default = \code{FALSE}.}

\item{ncore}{number of cores to run parallel computing Valid when \code{parallel = TRUE}.}

\item{...}{reserved passing parameter (or for future use).}
}
\value{
A data.frame containing mediation testing results of selected mediators.
\describe{
    \item{ID: }{Mediator ID/name.}
    \item{alpha: }{Coefficient estimates of exposure (X) --> mediators (M) (adjusted for covariates).}
    \item{beta: }{Coefficient estimates of mediators (M) --> outcome (Y) (adjusted for covariates and exposure).}
    \item{alpha*beta: }{The estimated indirect (mediation) effect of exposure on outcome through each mediator.}
    \item{rimp: }{Relative importance- the proportion of each mediator's mediation effect relative to the sum of the absolute mediation effects of all significant mediators.}
    \item{p-value: }{The joint p-value assessing the significance of each mediator's indirect effect, calculated based on the corresponding statistical approach.}
    \item{tau: }{The quantile level of the outcome (applicable only when using the quantile mediation model).}
}
}
\description{
\code{hima} is a wrapper function designed to perform various HIMA methods for estimating and testing high-dimensional mediation effects.
\code{hima} can automatically select the appropriate HIMA method based on the outcome and mediator data type.
}
\examples{
\dontrun{
# Note: In the following examples, M1, M2, and M3 are true mediators.

# Example 1 (continuous outcome - linear HIMA):
data(ContinuousOutcome)
pheno_data <- ContinuousOutcome$PhenoData
mediator_data <- ContinuousOutcome$Mediator

e1 <- hima(Outcome ~ Treatment + Sex + Age,
  data.pheno = pheno_data,
  data.M = mediator_data,
  mediator.type = "gaussian",
  penalty = "MCP", # Can be "DBlasso" for hima_dblasso
  scale = FALSE
) # Disabled only for simulation data
summary(e1)

# Efficient HIMA (only applicable to mediators and outcomes that are
# both continuous and normally distributed.)
e1e <- hima(Outcome ~ Treatment + Sex + Age,
  data.pheno = pheno_data,
  data.M = mediator_data,
  mediator.type = "gaussian",
  efficient = TRUE,
  penalty = "MCP", # Efficient HIMA does not support DBlasso
  scale = FALSE
) # Disabled only for simulation data
summary(e1e)

# Example 2 (binary outcome - logistic HIMA):
data(BinaryOutcome)
pheno_data <- BinaryOutcome$PhenoData
mediator_data <- BinaryOutcome$Mediator

e2 <- hima(Disease ~ Treatment + Sex + Age,
  data.pheno = pheno_data,
  data.M = mediator_data,
  mediator.type = "gaussian",
  penalty = "MCP",
  scale = FALSE
) # Disabled only for simulation data
summary(e2)

# Example 3 (time-to-event outcome - survival HIMA):
data(SurvivalData)
pheno_data <- SurvivalData$PhenoData
mediator_data <- SurvivalData$Mediator

e3 <- hima(Surv(Time, Status) ~ Treatment + Sex + Age,
  data.pheno = pheno_data,
  data.M = mediator_data,
  mediator.type = "gaussian",
  penalty = "DBlasso",
  scale = FALSE
) # Disabled only for simulation data
summary(e3)

# Example 4 (compositional data as mediator, e.g., microbiome):
data(MicrobiomeData)
pheno_data <- MicrobiomeData$PhenoData
mediator_data <- MicrobiomeData$Mediator

e4 <- hima(Outcome ~ Treatment + Sex + Age,
  data.pheno = pheno_data,
  data.M = mediator_data,
  mediator.type = "compositional",
  penalty = "DBlasso"
) # Scaling is always enabled internally for hima_microbiome
summary(e4)

#' # Example 5 (quantile mediation anlaysis - quantile HIMA):
data(QuantileData)
pheno_data <- QuantileData$PhenoData
mediator_data <- QuantileData$Mediator

# Note that the function will prompt input for quantile level.
e5 <- hima(Outcome ~ Treatment + Sex + Age,
  data.pheno = pheno_data,
  data.M = mediator_data,
  mediator.type = "gaussian",
  quantile = TRUE,
  penalty = "MCP", # Quantile HIMA does not support DBlasso
  scale = FALSE, # Disabled only for simulation data
  tau = c(0.3, 0.5, 0.7)
) # Specify multiple quantile level
summary(e5)
}

}
\references{
1. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J, Just A, Colicino E, Vokonas P, Zhao L,
Lv J, Baccarelli A, Hou L, Liu L. Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies.
Bioinformatics. 2016. DOI: 10.1093/bioinformatics/btw351. PMID: 27357171; PMCID: PMC5048064

2. Zhang H, Zheng Y, Hou L, Zheng C, Liu L. Mediation Analysis for Survival Data with High-Dimensional Mediators.
Bioinformatics. 2021. DOI: 10.1093/bioinformatics/btab564. PMID: 34343267; PMCID: PMC8570823

3. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L. Mediation Effect Selection in High-dimensional and Compositional Microbiome data.
Stat Med. 2021. DOI: 10.1002/sim.8808. PMID: 33205470; PMCID: PMC7855955

4. Zhang H, Chen J, Li Z, Liu L. Testing for Mediation Effect with Application to Human Microbiome Data.
Stat Biosci. 2021. DOI: 10.1007/s12561-019-09253-3. PMID: 34093887; PMCID: PMC8177450

5. Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L. HIMA2: High-dimensional Mediation Analysis and Its Application in
Epigenome-wide DNA Methylation Data. BMC Bioinformatics. 2022. DOI: 10.1186/s12859-022-04748-1. PMID: 35879655; PMCID: PMC9310002

6. Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L. High-Dimensional Quantile Mediation Analysis with Application to a Birth
Cohort Study of Mother–Newborn Pairs. Bioinformatics. 2024. DOI: 10.1093/bioinformatics/btae055. PMID: 38290773; PMCID: PMC10873903

7. Bai X, Zheng Y, Hou L, Zheng C, Liu L, Zhang H. An Efficient Testing Procedure for High-dimensional Mediators with FDR Control.
Statistics in Biosciences. 2024. DOI: 10.1007/s12561-024-09447-4.
}
