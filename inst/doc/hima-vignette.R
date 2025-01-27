## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE,
  echo = TRUE # Ensures all code is displayed by default
)
library(HIMA)

## ----hima-interface-----------------------------------------------------------
# hima(
#   formula,          # The model formula specifying outcome, exposure, and covariate(s)
#   data.pheno,       # Data frame with outcome, exposure, and covariate(s)
#   data.M,           # Data frame or matrix of high-dimensional mediators
#   mediator.type,    # Type of mediators: "gaussian", "negbin", or "compositional"
#   penalty = "DBlasso",  # Penalty method: "DBlasso", "MCP", "SCAD", or "lasso"
#   quantile = FALSE, # Use quantile mediation analysis (default: FALSE)
#   efficient = FALSE,# Use efficient mediation analysis (default: FALSE)
#   scale = TRUE,     # Scale data (default: TRUE)
#   sigcut = 0.05,    # Significance cutoff for mediator selection
#   contrast = NULL,  # Named list of contrasts for factor covariate(s)
#   subset = NULL,    # Optional subset of observations
#   verbose = FALSE   # Display progress messages (default: FALSE)
# )

## ----load-HIMA----------------------------------------------------------------
# library(HIMA)

## ----continuous-example-------------------------------------------------------
# data(ContinuousOutcome)
# hima_continuous.fit <- hima(
#   Outcome ~ Treatment + Sex + Age,
#   data.pheno = ContinuousOutcome$PhenoData,
#   data.M = ContinuousOutcome$Mediator,
#   mediator.type = "gaussian",
#   penalty = "MCP",
#   scale = FALSE # Demo data is already standardized
# )
# summary(hima_continuous.fit, desc=TRUE)
# # `desc = TRUE` option to show the description of the output results

## ----efficient-example--------------------------------------------------------
# hima_efficient.fit <- hima(
#   Outcome ~ Treatment + Sex + Age,
#   data.pheno = ContinuousOutcome$PhenoData,
#   data.M = ContinuousOutcome$Mediator,
#   mediator.type = "gaussian",
#   efficient = TRUE,
#   penalty = "lasso",
#   scale = FALSE # Demo data is already standardized
# )
# summary(hima_efficient.fit, desc=TRUE)
# # Note that the efficient HIMA is controlling FDR

## ----binary-example-----------------------------------------------------------
# data(BinaryOutcome)
# hima_binary.fit <- hima(
#   Disease ~ Treatment + Sex + Age,
#   data.pheno = BinaryOutcome$PhenoData,
#   data.M = BinaryOutcome$Mediator,
#   mediator.type = "gaussian",
#   penalty = "MCP",
#   scale = FALSE # Demo data is already standardized
# )
# summary(hima_binary.fit)

## ----survival-example---------------------------------------------------------
# data(SurvivalData)
# hima_survival.fit <- hima(
#   Surv(Time, Status) ~ Treatment + Sex + Age,
#   data.pheno = SurvivalData$PhenoData,
#   data.M = SurvivalData$Mediator,
#   mediator.type = "gaussian",
#   penalty = "DBlasso",
#   scale = FALSE # Demo data is already standardized
# )
# summary(hima_survival.fit)

## ----microbiome-example-------------------------------------------------------
# data(MicrobiomeData)
# hima_microbiome.fit <- hima(
#   Outcome ~ Treatment + Sex + Age,
#   data.pheno = MicrobiomeData$PhenoData,
#   data.M = MicrobiomeData$Mediator,
#   mediator.type = "compositional",
#   penalty = "DBlasso"
# )
# summary(hima_microbiome.fit)

## ----quantile-example---------------------------------------------------------
# data(QuantileData)
# hima_quantile.fit <- hima(
#   Outcome ~ Treatment + Sex + Age,
#   data.pheno = QuantileData$PhenoData,
#   data.M = QuantileData$Mediator,
#   mediator.type = "gaussian",
#   quantile = TRUE,
#   penalty = "MCP",
#   tau = c(0.3, 0.5, 0.7),
#   scale = FALSE # Demo data is already standardized
# )
# summary(hima_quantile.fit)

