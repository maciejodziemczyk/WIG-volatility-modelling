# WIG volatility modelling
Project created for *Financial Markets Modelling* (org. Modelowanie Rynków Finansowych (IiE)) classes at WNE UW

Language:
 - Polish - classes, report and code comments

Semester: III (MA studies)

## About
The main objective of this project was to do own research based on literature (given paper). Bartek Kuźma and I modelled WIG volatility on daily and weekly data. Our research was about GARCH and naive models comparison. We chose several GARCH family models s.t standard GARCH, Exponential Garch (E-GARCH), Threshold Garch (T-GARCH) and Component GARCH (C-GARCH) with different epsilon distribution assumptions (Normal, t-Student, skewed t-Student and Generalized Error Distribution). Considered naive model were Random Walk, Historical Average and Moving Average (all coded from scratch). We compared our modes using many metrics s.t Mean Error, Mean Absolute Error, Root Mean Square Error, Adjusted Mean Absolute Percentage Error, Thein Income Coefficient (symmetric metrics), MME(U) and MME(O) (asymmetric metrics), DCP and DCP(U) (financial loss functions).

We expected more sophisticated methods to be better that the simple ones, moreover we expected GARCHs with asymmetric epsilon distributions to be better that with the symmetrics ones due to financial instruments distributions properties (leptokurthosis and skewness). 

In this research we downloaded WIG data from stooq.pl and performed standard preprocessing steps - log of the first differences. We also performed standard data analysis s.t. ACF and PACF functions, ARCH LM tests, Durbin-Watson test, empirical distribution analysis (histogram and kernel distribution esitmation) and basic dsscriptive statistics.

To compare methods we simulated stochastic process with for loop and models reestimation with rolling window, so we obtained out of sample results (whole 2020).

Findings
 - there was no difference between GARCH models with different distribution assumptions 
 - GARCH models was better on higher frequency data (daily) 
 - for weekly data GARCH models were more likely to overestimate while naive models were more likely to underestimate volatility
 - results variance was lower for higher frequency data
 - no one the best model, different models won depend on metric

## Repository Description
 - Data-processing.ipynb - preprocessing with python (daily data) - we found that we could do this easily in R so weekly data was preprocessed in R :)
 - GarchOOSresults.csv - models results for daily data
 - Kuzma_Odziemczyk_model_prezentacja.pdf - our project presentation (we had to present it to the group and teacher)
 - Kuzma_Odziemczyk_model_raport.pdf - our project report (our asignment), I encourage you to read it as always ;)
 - MRF.R - R script with analysis (easy to deal with via R studio and R project enviroment)
 - other folders .csv and txt are just a datasets and our dumped results

## Technologies
 - R (analysis)
 - LaTeX (presentation and report)

## Authors
 - Bartłomiej Kuźma
 - Maciej Odziemczyk
