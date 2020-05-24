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

To access the AugBin GUI, go to https://martinamcm.shinyapps.io/multsampsize/. 

An example dataset for a composite with one continuous and one binary component is available in the repository.

## Tutorial

### Co-Primary and Multiple Primary Endpoints



