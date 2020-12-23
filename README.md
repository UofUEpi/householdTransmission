The R code in these files produce the results in the manuscript "High variability in transmission of SARS-CoV-2 within households and implications for control" 

"analysisData.txt" contains the household data used by the analyses in the manuscript. The column definitions are described in the main text of the manuscript.

"Table1MLE.R" contains R code that produces the main maximum likelihood results in Table 1 of the manuscript. It calls the code in "prepareData.R" and "MLEfunctions.R" that calculate the likelihood against the data. These codes are also called within other R codes described below.

"Table1Bootstrap.R” contains R code for the bootstrap estimates in Table 1 of the manuscript. It makes use of “bootstrapResults.txt” that contains MLE fits to each of the 500 simulated data sets contained in the zipped folder “SimData.zip.” The code that calculated those MLE for the simulated data sets is found in “fitMLEtoSimData.R” and the code that stochastically simulates synthetic data sets is found in “generateSimData.R”
 
“Table2.R” and “Table3.R” contain R code that produces the results in Tables 2 and 3 of the manuscript.

“gammaDispersionResults.R” contains R code that calculates the result described in the manuscript for the dispersion parameter k (shape parameter of gamma distribution) that matches the mean and variance of the MLE household transmission probability.

“householdRresults.R” contains R code that calculates the household reproduction numbers Rh and Rh* and threshold values of non-household reproduction number Rc, as described in the manuscript results.

“SupplementFigures.R” contains R code that produces all the figures in the supplementary results of the manuscript. Some plotting display commands may be specific to R for Windows.

 “TableS1.R” and “TableS2.R” contain R code that produces the results in Tables S1 and S2 in the supplementary material of the manuscript.


