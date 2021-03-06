% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelSpecificationAndEstimation.R
\name{modLik}
\alias{modLik}
\alias{model_translateParameters}
\alias{model_makeDefaultParameterStructures}
\alias{model_fellerConditionCheck}
\alias{model_Likelihood_extraNoise}
\alias{model_wrapLikelihood_extraNoise}
\alias{model_Likelihood_portfolio_extraNoise}
\alias{model_wrapLikelihood_portfolio}
\alias{model_Likelihood_affineContract}
\alias{model_wrapLikelihood_affineContract}
\title{Affine model likelihood based on divergence prices}
\usage{
model_translateParameters(par.vec, par.names = names(par.vec), par.restr,
  N.factors)

model_makeDefaultParameterStructures(N.factors,
  pq.equality = c("Q$jmp$lvec", paste0("Q$jmp$lprop.", 1:N.factors),
  paste0("Q$", 1:N.factors, "$eta")))

model_fellerConditionCheck(params.P, params.Q, N.factors)

model_Likelihood_extraNoise(data.structure, model.spec,
  for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5,
  penalized = FALSE, penalty = 1e+12)

model_wrapLikelihood_extraNoise(data.structure, model.spec,
  for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5,
  penalized = FALSE, penalty)

model_Likelihood_portfolio_extraNoise(data.structure, model.spec,
  for.estimation = FALSE,
  filterFoo = divergenceModelR:::portfolio_sqrtFilter, N.points = 5,
  penalized = FALSE, penalty = 1e+12, N.GL.points = 96)

model_wrapLikelihood_portfolio(data.structure, model.spec,
  for.estimation = FALSE,
  filterFoo = divergenceModelR:::portfolio_sqrtFilter, N.points = 5,
  penalized = FALSE, penalty, N.GL.points = 96)

model_Likelihood_affineContract(data.structure, model.spec,
  for.estimation = FALSE, filterFoo = DSQ_sqrtFilter, N.points = 5,
  penalized = FALSE, penalty = 1e+12)

model_wrapLikelihood_affineContract(data.structure, model.spec,
  for.estimation = FALSE,
  filterFoo = divergenceModelR:::portfolio_sqrtFilter, N.points = 5,
  penalized = FALSE, penalty)
}
\arguments{
\item{par.vec}{vector with model parameter values}

\item{par.names}{parameter names, character vector equal in length to par}

\item{par.restr}{parameter equality restrictions, data.frame; par.vec and par.restr have to exhaust the model parameter set together.}

\item{N.factors}{integer, number of SV factors}

\item{data.structure}{\code{list} with fields \code{spec.mat}, \code{obs.data}, \code{noise.cov.cube}}

\item{model.spec}{\code{list} with fields \code{params.P}, \code{params.Q}, \code{jump.type}, \code{dt}, \code{N.factors}, \code{error}, \code{mkt}}

\item{for.estimation}{\code{logical}; determines return type (log-lik) or filtering result}

\item{filterFoo}{\code{function} that handles the filtering, must correspond to model specification and to provided observables.}

\item{N.points}{\code{integer}, number of integration points for double quadrature in moments of the state and stock price.}

\item{penalized}{\code{FALSE} by default, if \code{TRUE}, the Feller constraint is imposed as a penalty on the likelihood. Otherwise, infinity is returned (discontinuous at boundary).}

\item{noisePar}{vector of noise variance magnitudes, equal to number of observed pfolios}

\item{noisePar}{vector of noise variance magnitudes, equal to number of observed pfolios}

\item{noisePar}{vector of noise variance magnitudes, equal to number of observed pfolios}
}
\value{
\code{model_Likelihood} if \code{for.estimation==TRUE}: log-likelihood value (NOT negative of...), else: list with filtering results

\code{model_translateParameters} \code{list} with fields \code{P} and \code{Q}, input for all ODE calling functions.

\code{model_makeDefaultParameterStructures} returns \code{data.frame} \code{par.restr} and \code{character} vector \code{par.names}

\code{model_fellerConditionCheck} list with two logical vectors reporting whether the Feller conditions are satisfied

\code{model_Likelihood_extraNoise} \code{list} with fields \code{P} and \code{Q}, input for all ODE calling functions.

\code{model_wrapLikelihood_extraNoise} wraps the likelihood function with extra noise so that it only accepts a parameter vector argument -- use this for optimizers that do not allow passing extra arguments to the optimised function.

\code{model_Likelihood_portfolio_extraNoise} \code{list} with fields \code{P} and \code{Q}, input for all ODE calling functions.

\code{model_wrapLikelihood_portfolio} wraps the likelihood function so that it only accepts a parameter vector argument -- use this for optimizers that do not allow passing extra arguments to the optimised function.

\code{model_Likelihood_affineContracts} \code{list} with fields \code{P} and \code{Q}, input for all ODE calling functions.

\code{model_wrapLikelihood_affineContract} wraps the likelihood function so that it only accepts a parameter vector argument -- use this for optimizers that do not allow passing extra arguments to the optimised function.
}
\description{
Functions for preparing model specification, evaluating the likelihood for the DSQ filtering.
}
\details{
Not much for now
}
