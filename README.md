# Dynamic Estimation and Forecasting of the Euro Area Yield Curve

This repository provides the full set of MATLAB scripts used for the dynamic estimation and forecasting of the Euro Area yield curve using the Nelson-Siegel model framework. The methodology developed and implemented here is part of a **master's thesis in Quantitative Finance** and builds on seminal works by Diebold and Li (2006) as well as Diebold, Rudebusch, and Aruoba (2006).


## 🔍 Overview
The dataset used spans from **September 2004 to January 2025**, and is divided into two subsamples:
- **Training set**: September 2004 – February 2017  
- **Validation set**: March 2017 – January 2025


## 🛠️ Models Implemented for Estimation
The Nelson-Siegel model is applied in two distinct estimation settings:

1. **State-space representation with Kalman filter**:  
   The model is reformulated in state-space form, where the latent factors (level, slope, and curvature) evolve according to a VAR(1) process. Estimation of latent states is carried out using a Kalman filter for real-time estimates and a Kalman smoother for full-sample inference. The model also incorporates parameter estimation via maximum likelihood.

2. **Two-step estimation via autoregressive models**:  
   Latent factors are first estimated by minimizing the distance between observed and model-implied yields at each time point. Subsequently, an autoregressive process is fit to the estimated time series of factors. This variant follows the traditional Diebold-Li (2006) estimation strategy.


## 📊 Forecasting Methodology
The core component of this project lies not only in the estimation of yield curves, but also in the recursive **out-of-sample forecasting** of future yield curves. Forecasting is performed using an expanding window procedure. At each point in the validation period, models are re-estimated using all available data up to that date, and forecasts are generated for multiple horizons. This approach is consistent with the methods used in the above-referenced academic literature.

In addition to the two Nelson-Siegel variants, the code includes benchmark forecasting models for comparison:
- AR(1) processes applied directly to yields at each maturity
- VAR(1) applied to the vector of yields
- A random walk (no-change) model, which serves as the benchmark in forecast evaluation

To assess and compare the forecast accuracy of the models, the **Diebold-Mariano test** is implemented. This statistical test allows for pairwise comparison of predictive performance across models and maturities, taking into account forecast errors over time. In particular, it evaluates whether the predictive accuracy of a model is significantly different from that of the random walk or alternative competing models.


## 📁 What's Inside
- `/DNS.m/`: Main and function codes, provided in Matlab
- `/ECBData.xlsx/`: Dataset provided by the European Central Bank (ECB), converted to a monthly frequency for the analysis conducted in this thesis (3-month to 10-year maturities).
- `/AR_factor_forecasting/`: For further details regarding the specifications of the models used for comparison with the random walk, please contact the author via his LinkedIn profile (_**below in Credits**_). The process is quite the same for all of them.
- `/Diebold_Mariano_test/`: Implementation of Diebold-Mariano test across models and maturities.


## 🚀 Requirements
- MATLAB R2024a or later
- Econometrics Toolbox
- Optimization Toolbox


## 🎓 Credits
For further information, clarifications, or inquiries regarding the methodology, dataset, or code, feel free to reach out to the repository maintainer. Additional details can be found on the LinkedIn profile:

[Davide Di Virgilio's LinkedIn Profile](https://www.linkedin.com/in/davide-di-virgilio-b904a2261)


## 🧾 Citation

This repository is part of a **Master’s Thesis in Quantitative Finance**.   
For academic use, please cite the original papers:

> Diebold, F.X., & Li, C. (2006). Forecasting the Term Structure of Government Bond Yields. *Journal of Econometrics*.  
> Diebold, F.X., Rudebusch, G.D., & Aruoba, S.B. (2006). The Macroeconomy and the Yield Curve: A Dynamic Latent Factor Approach.


 


