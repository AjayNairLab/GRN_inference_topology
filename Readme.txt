This repository gives data and codes used for generating results for the Journal paper:
*******************************************************************
Improving gene regulatory network inference using network topology information; A. Nair, M. Chetty, and P. P. Wangikar, Mol. BioSyst., 2015, DOI: 10.1039/C5MB00122F.
*******************************************************************

If you find these useful in your work, please cite the above paper.

Some of the codes are from the toolbox called GlobalMIT. These code files have that information explicitly mentioned. Some of the codes are modified from GlobalMIT and these are also mentioned at the beginning of the files. Others have been developed me. I am making this available as a requirement for the above journal review process; and also because the results in the paper show a good improvement in the dynamic Bayesian network inference computational time. Thus, these codes will be a good supplement for the details of the algorithm given in the paper.

(c) 2014-2015 Ajay Nair

The details of the files present in this repository is given below:
Algorithms:
1. Normal: globalMIT_ab.m
2. maxP: globalMIT_ab_maxP.m
3. maxPiter: globalMIT_ab_maxP_iter2.m
4. maxPincrement: globalMIT_ab_maxP_incr.m

Other required functions:
conditional_MI_DBN_ab.m : for calculating the conditional MI
findLexicalIndex.m: finding the lexical index for saving intermediate results of inference scores
fnPerformanceMeasure.m: to caculate the different parameters of inference performance like True positives, true negatives, false positives, false negatives, precision, recall, f-score, and specificity
multi_time_series_cat.m: to concatenate multiple time series data.
myIntervalDiscretize.m: to discretize continious data

Tutorial files: E. coli SOS network.
data_samples_SOS.mat: the data file for E. coli SOS network data
net_SOS_normal.m: normal network inference with globalMIT
net_SOS_GMIT_maxP.m: network inference of SOS network with maxP algorithm
net_SOS_maxPiter.m: network inference of SOS network with maxPiter algorithm
net_SOS_maxPincrement.m: network inference of SOS network with maxPincrement algorithm
