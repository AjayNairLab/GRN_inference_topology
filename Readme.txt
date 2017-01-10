This repository gives data and codes used for generating results for the Journal paper:
*******************************************************************
Nair, A., Chetty, M., and Wangikar, P.P. (2015). Improving gene regulatory network inference using network topology information. Mol. BioSyst. 11, 2449–2463.
DOI: 10.1039/C5MB00122F.
*******************************************************************

If you find these useful in your work, please cite the above paper.

Some of the codes are from the toolbox called GlobalMIT (http://code.google.com/p/globalmit/). These code files have that information explicitly mentioned. Some of the codes are modified from GlobalMIT and these are also mentioned at the beginning of the files. Others have been developed by me. I am making this available as a requirement for the above journal review process; and also because the results in the paper show a good improvement in the dynamic Bayesian network inference computational time. Thus, these codes will be a good supplement to the details of the algorithm given in the paper. Further, the results can be easily verified and reproduced.

(c) 2014-2015 Ajay Nair

The details of the files present in this repository are given below:

Algorithms:
-----------
1. Normal: globalMIT_ab.m
2. maxP: globalMIT_ab_maxP.m
3. maxPiter: first iteration performed with 'globalMIT_ab_maxP.m' and the second iteration performed with 'globalMIT_ab_maxP_iter2.m'
4. maxPincrement: first iteration performed with 'globalMIT_ab_maxP.m' and the second and subsequent iterations are performed with 'globalMIT_ab_maxP_incr.m'


Other required functions:
-------------------------
conditional_MI_DBN_ab.m : for calculating the conditional MI
findLexicalIndex.m: finding the lexical index for saving intermediate results of inference scores
fnPerformanceMeasure.m: to caculate the different parameters of inference performance like True positives, true negatives, false positives, false negatives, precision, recall, f-score, and specificity
multi_time_series_cat.m: to concatenate multiple time series data.
myIntervalDiscretize.m: to discretize continious data


Tutorial files: E. coli SOS network.
------------------------------------
data_samples_SOS.mat: the gene expression data file for E. coli SOS network from Ronen et al., 2003. It contains 2 variables:
  1. sosPromoterAct: 4x8x50 matrix with 4-experiments; 8-genes; and 50-samples
  2. time: 1x5o matrix. These are the time intervals of 50 samples  
net_SOS_normal.m: normal network inference with globalMIT
net_SOS_maxP.m: network inference of SOS network with maxP algorithm
net_SOS_maxPiter.m: network inference of SOS network with maxPiter algorithm
net_SOS_maxPincrement.m: network inference of SOS network with maxPincrement algorithm


Data files for the other networks.
----------------------------------
IRMA:
	data_samples_irma_SOffDct.mat; %Switch-off time series data set
	data_samples_irma_SOnDct.mat; %Switch-on time series data set
	net_irma_normal.m; %normal network inference with globalMIT
	net_irma_maxP.m; %network inference with maxP algorithm
	net_irma_maxPiter.m; %network inference with maxPiter algorithm
	net_irma_maxPincrement.m; %network inference with maxPincrement algorithm
	
DMn10e11:
	data_samples_DMn10e11.mat; %gene expression data file
	net_DMn10e11_maxPincrement.m; %network inference maxPincrement algorithm
DMn10e15
	data_samples_DMn10e15.mat; %gene expression data file
	net_DMn10e15_maxPincrement.m; %network inference with maxPincrement algorithm
GNWn20e50
	data_samples_GNWn20e50.mat; %gene expression data file
	net_GNWn20e50_maxPincrement.m; %network inference with maxPincrement algorithm
GNWn20e33
	data_samples_GNWn20e33.mat; %gene expression data file
	net_GNWn20e33_maxPincrement.m; %network inference with maxPincrement algorithm
GNWn20e31
	data_samples_GNWn20e31.mat; %gene expression data file
	net_GNWn20e31_maxPincrement.m; %network inference with maxPincrement algorithm
GNWn20e19
	data_samples_GNWn20e19.mat; %gene expression data file
	net_GNWn20e19_maxPincrement.m; %network inference with maxPincrement algorithm
n10e9
	data_samples_n10e9.mat; %gene expression data file
	net_n10e9_maxPincrement.m; %network inference with maxPincrement algorithm
n10e26
	data_samples_n10e26.mat; %gene expression data file
	net_n10e26_maxPincrement.m; %network inference with maxPincrement algorithm
n10e45
	data_samples_n10e45.mat; %gene expression data file
	net_n10e45_maxPincrement.m; %network inference with maxPincrement algorithm
n10e9LS
	data_samples_n10e9.mat; %gene expression data file
	net_n10e9LS_maxPincrement.m; %network inference with maxPincrement algorithm
n10e26LS
	data_samples_n10e26.mat; %gene expression data file
	net_n10e26LS_maxPincrement.m; %network inference with maxPincrement algorithm
n10e45LS
	data_samples_n10e45.mat; %gene expression data file
	net_n10e45LS_maxPincrement.m; %network inference with maxPincrement algorithm
YUn20e9
	data_samples_YUn20e9.mat; %gene expression data file
	net_YUn20e9_maxPincrement.m; %network inference with maxPincrement algorithm

n20e19
	data_samples_n20e19.mat; %gene expression data file
	net_n20e19_maxPincrement.m; %network inference with maxPincrement algorithm
n20e36
	data_samples_n20e36.mat; %gene expression data file
	net_n20e36_maxPincrement.m; %network inference with maxPincrement algorithm
n20e150
	data_samples_n20e150.mat; %gene expression data file
	net_n20e150_maxPincrement.m; %network inference with maxPincrement algorithm

