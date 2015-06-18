These are the codes used for generating the results in the paper.

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
