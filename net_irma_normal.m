%This file is part of the data and codes used for generating results for the Journal paper:
%*******************************************************************
%Improving gene regulatory network inference using network topology information; 
%A. Nair, M. Chetty, and P. P. Wangikar, Mol. BioSyst., 2015, DOI: 10.1039/C5MB00122F.
%*******************************************************************

%If you find these useful in your work, please cite the above paper.
%(c) 2014-2015 Ajay Nair

%This file is a tutorial on IRMA network reconstruction using normal inference 
%with globalMIT algorithm (a optimal DBN inference algorithm)

%Usage:
% net_irma_normal()

function []=net_irma_normal()
%*********************CONSTANTS REQUIRED
n_state=3;%no of discrete states for the microarray data
alpha=0.999;%significance level for the mutual information test for independance.
allowSelfLoop=0;%allow self regulated link (=1) or not (=0)
%**************************************
%clc
%**************************************Input data
% switch ON data
load data_samples_irma_SOnDct.mat; % load the time series data of 2^-(DeltaCt)
%data should be in the format [experiments X genes]
a1d=myIntervalDiscretize(a1,n_state); %descretize according to the rows
a2d=myIntervalDiscretize(a2,n_state); %descretize according to the rows
a3d=myIntervalDiscretize(a3,n_state); %descretize according to the rows
a4d=myIntervalDiscretize(a4,n_state); %descretize according to the rows
a5d=myIntervalDiscretize(a5,n_state); %descretize according to the rows

%switch OFF data
load data_samples_irma_SOffDct.mat; % load the different time series data in IRMA paper supplement
%data should be in the format [experiments X genes]
b1d=myIntervalDiscretize(b1,n_state); %descretize according to the rows
b2d=myIntervalDiscretize(b2,n_state); %descretize according to the rows
b3d=myIntervalDiscretize(b3,n_state); %descretize according to the rows
b4d=myIntervalDiscretize(b4,n_state); %descretize according to the rows
[b,c]=multi_time_series_cat(a1d,a2d,a3d,a4d,a5d,b1d,b2d,b3d,b4d);% combining both S-ON and S-OFF
%network nodes
%SWI5 = 1; CBF1 = 2; GAL4 = 3; GAL80 = 4; ASH1=5;
nodeNames=[{'SWI5'},{'CBF1'},{'GAL4'},{'GAL80'},{'ASH1'}];
%the actual IRMA network         
actualNet=[ 0  1  0  1  1;
            0  0  1  0  0;
            1  0  0  0  0;
            0  0  0  0  0;
            0  1  0  0  0];
%**************************************

%**************************************Performing inference
 tic();
  [best_net]=globalMIT_ab(b,c,alpha,allowSelfLoop);
 t=toc();
 
fprintf('\nActual Network:\n');
actualNet

fprintf('Learned Network:\n');
best_net

fprintf('Performance Measures:\n True positives, true negatives, false positives, false negatives, precision, recall, f-score, specificity \n');
%tp tn fp fn prec recl fscor spec
%DOING PERFORMANCE MEASURES
M=fnPerformanceMeasure(best_net, actualNet)

fprintf('Time taken in seconds: %f\n',t);



%creating the network figure using 'createDotGraphic' function
%createDotGraphic(actualNet,nodeNames,'Original SOS network');
%createDotGraphic(best_net,nodeNames,'Learned SOS network');


end
