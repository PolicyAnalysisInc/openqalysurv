# openqalysurv

[![CircleCI build
status](https://circleci.com/gh/PolicyAnalysisInc/openqalysurv.svg?style=svg)](https://circleci.com/gh/PolicyAnalysisInc/openqalysurv)
[![Codecov test
coverage](https://codecov.io/gh/PolicyAnalysisInc/openqalysurv/branch/master/graph/badge.svg)](https://codecov.io/gh/PolicyAnalysisInc/openqalysurv?branch=master)

## Overview

openqalysurv is a grammar of survival modeling for use with openqaly and
the openqalyverse. It provides a series of verbs that help you define,
combine, and modify survival distributions.

## Installation

``` r

# Install latest version from github
remotes::install_github("PolicyAnalysisInc/openqalysurv")
```

## Defining Survival Distributions

- Parametric distributions with specified parameter values:
  [`define_surv_param()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_param.md)
- Parametric mixture & non-mixture cure models with specified parameter
  values:
  [`define_surv_cure()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_cure.md)
- Royston & Parmar spline models with specified parameter values:
  [`define_surv_spline()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_spline.md)
- Life-tables:
  [`define_surv_lifetable()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_lifetable.md)
- Kaplan-Meiers based on product-limit table:
  [`define_surv_km()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_km.md)
- Custom Functions:
  [`define_surv_func()`](https://policyanalysisinc.github.io/openqalysurv/reference/define_surv_func.md)
- Full support for models estimated using
  [flexsurv](https://cran.r-project.org/web/packages/flexsurv/index.html)
  package
- Partial support for KMs & Cox Models estimated using
  [survival](https://cran.r-project.org/web/packages/survival/index.html)
  package

## Modifying Survival Distributions

- Apply hazard ratio:
  [`apply_hr()`](https://policyanalysisinc.github.io/openqalysurv/reference/apply_hr.md)
- Apply odds ratio:
  [`apply_or()`](https://policyanalysisinc.github.io/openqalysurv/reference/apply_or.md)
- Apply acceleration factor:
  [`apply_af()`](https://policyanalysisinc.github.io/openqalysurv/reference/apply_af.md)
- Apply shift in time:
  [`apply_shift()`](https://policyanalysisinc.github.io/openqalysurv/reference/apply_shift.md)
- Set covariate levels of models with covariates:
  [`set_covariates()`](https://policyanalysisinc.github.io/openqalysurv/reference/set_covariates.md)

## Combining Survival Distributions

- Join survival distributions together at specified times:
  [`join()`](https://policyanalysisinc.github.io/openqalysurv/reference/join.md)
- Mix together survival distributions with specified weights:
  [`mix()`](https://policyanalysisinc.github.io/openqalysurv/reference/mix.md)
- Combine distributions as independent risks:
  [`add_hazards()`](https://policyanalysisinc.github.io/openqalysurv/reference/add_hazards.md)

## Generating Predicted Values

- Generate survival probabilities:
  [`surv_prob()`](https://policyanalysisinc.github.io/openqalysurv/reference/surv_prob.md)
- Generate event probabilities:
  [`event_prob()`](https://policyanalysisinc.github.io/openqalysurv/reference/event_prob.md)
