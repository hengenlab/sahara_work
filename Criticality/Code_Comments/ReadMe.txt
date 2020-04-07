Check AvalancheAnalysisExample.m. This script shows an example on a dataset from deprived hemisphere during baseline period.

Data includes the dataset that already be transferred to binary data.
each column is a time bin, while each row is a unit.

AV_analysis_BurstT.m is used to extract avalanche sizes and durations. The only parameter need to tune is 'perc', which determines where could we set the threshold.

AV_analysis_ExponentErrorComments.m is applied to get the power fitting. Please check the comments in AvalancheAnalysisExample.m for details.

tplfit.m is used to fit a truncated power law with log maximum likelihood method. plmle.m is also for exponents fitting.


EXCLUDE.m is used to extract the lower and upper boundaries.

pvaluenew.m generate Niter(1000) surrogated datasets and compared with experimental data to return us a pvalue. When pvalue is large (>0.05), the null hypothesis could not be rejected and the experimental follows power law good enough.