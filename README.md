# MultSampSize

## [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)

## Description
MultSampSize is a Shiny app which can be used to calculate the sample size required for studies using a latent variable to analyse mixed outcome endpoints in the form of:

* Co-primary - requires an effect in all outcomes
* Multiple primary - requires an effect in at least one outcome
* Composite - requires an effect in some overall combination of the outcomes

The latent variable model employed for the mixed endpoints assumes the discrete variables are manifestations of latent continuous variables. The observed continuous and latent continuous outcomes are assumed to follow a multivariate normal distribution. More details on sample size determination using these models is available at <https://arxiv.org/abs/1912.05258>.

The tutorial below provides step-by-step guidance for using the MultSampSize Shiny app. In the case that further queries arise about the functionality of the app for specific applications, contact Martina McMenamin at <martina.mcmenamin@mrc-bsu.cam.ac.uk>.

## Getting started

To access the MultSampSize GUI, go to https://martinamcm.shinyapps.io/multsampsize/. 

An example dataset for a composite with one continuous and one binary component is available in the repository.

## Tutorial

The underlying model assumed for analysis is the same for each of the endpoints considered, however it is employed in different ways. The left side bar on the landing page allows the user to choose the relevant endpoint from co-primary, multiple primary and composite endpoints as shown. 

<p align="center">
<img src="/Images/SampleSizeHome.png" title="SampleSizeHome" width="100%" />
</p>

### Co-Primary and Multiple Primary Endpoints

The co-primary case will be demonstrated through an example, however the functionality also applies to the multiple outcome scenario. The user should begin by choosing the number of continuous and binary outcomes that make up the endpoint, which can be one or two continuous and zero or one binary measures. Clicking the 'Generate Model' button shows a summary of the latent variable model used and the power function calculated. This is shown below for the case when the co-primary endpoint is made up of two continuous and one binary outcome.

<p align="center">
<img src="/Images/CoPrim1.png" title="CoPrim1" width="100%" />
</p>

Note that the model employed is the same for the co-primary and multiple primary endpoints however the power function differs.

<p align="center">
<img src="/Images/MultPrim1.png" title="MultPrim1" width="100%" />
</p>

When the structure of the endpoint has been selected the user can input parameter values where &delta;<sub>k</sub> is the risk difference, &pi;<sub>Tk</sub> and &pi;<sub>Ck</sub> is the probability that observed binary outcome Y<sub>ik</sub> is equal to  1 for patient i, &sigma;<sub>k</sub> is the standard deviation in outcome k and &rho;<sub>kl</sub>} is the correlation between outcome k and l. 

<p align="center">
<img src="/Images/CoPrim2.png" title="CoPrim2" width="100%" />
</p>

The alpha level, desired power and maximum number of subjects can be adjusted in the 'Sample Size Estimation' panel as shown below. The resulting power curve is displayed and is updated when these or the model parameters are adjusted. The number of subjects required per arm is stated below the plot, as shown. 

<p align="center">
<img src="/Images/CoPrim3.png" title="SampleSize" width="100%" />
</p>


### Composite Endpoints

For composite responder endpoints, the user must select the number of continuous and binary components along with the dichotomisation threshold for each of the continuous measures. 

<p align="center">
<img src="/Images/Comp1.png" title="Comp1" width="40%" />
</p>

After selecting 'Get Model' a model summary, treatment effect definition and the power function will be displayed in the 'Model Summary' panel.

<p align="center">
<img src="/Images/Comp2.png" title="Comp2" width="90%" />
</p>

The user must upload pilot data and click 'Obtain estimates' to provide parameter estimates using the panel as shown. A loading bar will be shown when the model is being fitted and a table with parameter estimates from both the latent approach and standard binary approach will be displayed when it is complete. 

<p align="center">
<img src="/Images/Comp3.png" title="Comp3" width="80%" /><img src="/Images/Comp4.png" title="Comp4" width="80%" />
</p>

The uploaded data must have columns ordered as follows:
* Patient id
* Observed continuous outcomes 
* Observed binary outcome (if any)
* Baseline measures for observed continuous measures 

The 'Sample Size Estimation' panel shows a plot of the power curve for both the latent variable model and the standard binary model, using the &delta; and &sigma; estimates from the panel above. The one-sided significance level, power target and number of patients on the x-axis can be altered using the controls shown.

<p align="center">
<img src="/Images/Comp5.png" title="Comp5" width="100%" />
</p>
