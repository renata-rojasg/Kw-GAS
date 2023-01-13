<div>
  <h1> The Kumaraswamy generalized autoregressive model for the quantiles of double-bounded hydro-environmental time series </h1> 
</div>

This repository contains an R implementation of the Kumaraswamy generalized autoregressive (Kw-GAS) model. 

<h2> Available implementations and files: </h2>

<h3> To fit the models: </h3>

<ul>
  <li> kum-q-phi.R: auxiliary file with the probability density function (dkum), 
  cumulative density function (pkum), 
  quantile function (qkum), 
  and a function for generating pseudo-random numbers (rkum) of the Kumaraswamy distribution with the quantile-based parametrization.</li>
  <li>GASKwFit-COV.R: function GASKw.fit() estimates the Kw-GAS.</li>
  <li>GASBetaFit-COV.R: function GASbeta.fit() estimates the beta-GAS.</li>  
</ul>  

<h3> To simulate the models: </h3>

<ul>
  <li>simu_Beta-GAS-COV.R: function simu.GASKw() generates a sample of the Kw-GAS model.</li>
  <li>simu_Beta-GAS-COV.R: function simu.GASbeta() generates a sample of the beta-GAS model.</li>  
  <li> simu_Table1.R: simulation study to evaluate the parameter estimatiors of the Kw-GAS model.</li>
  <li> simu_Table2.R: simulation study to evaluate the diagnostic measures.</li>
</ul> 

<h3> Water reservoirs analysis: </h3>

<ul>
  <li> reservoirs_data.csv: data set .</li>
  <li> app_reservatorios-FITTED.R: implementations for the fitted models and related measures.</li>  
  <li> app_reservatorios-FORECAST.R: implementations for out-of-sample forecasts and related measures.</li>  
  <li> app_reservatorios.RData: results for the fitted models and related measures.</li>  
  <li> app_reservatorios-FORECAST1step.RData: results for out-of-sample forecasts and related measures.</li>  
  <li> dashboard_resevoirs.Rmd: codes for interactive dashboard is available at https://renatarguerra.shinyapps.io/UV-Reservoirs/.</li>
</ul> 


Created by Renata Rojas Guerra (renata.r.guerra@ufsm.br)
