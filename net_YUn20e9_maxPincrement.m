%This file is part of the data and codes used for generating results for the Journal paper:
%*******************************************************************
%Improving gene regulatory network inference using network topology information; 
%A. Nair, M. Chetty, and P. P. Wangikar; Mol. BioSyst., 2015, DOI: 10.1039/C5MB00122F.
%*******************************************************************

%If you find these useful in your work, please cite the above paper.
%(c) 2014-2015 Ajay Nair

%This file is a tutorial on YUn20e9 network reconstruction using maxPincrement algorithm

%Usage:
% net_YUn20e9_maxPincrement()
function []=net_YUn20e9_maxPincrement()
%*********************CONSTANTS REQUIRED
n_state=3;%no of discrete states for the microarray data
alpha=0.999;%significance level for the mutual information test for independance.
allowSelfLoop=0;%allow self regulated link (=1) or not (=0)
maxParent=1; %maximum number of parents a node can have
maxParentLimit=10;
%**************************************
%clc
%**************************************Input data
load data_samples_YUn20e9 data; % load the time series data 
%data should be in the format [exps X genes]
a1d=myIntervalDiscretize(data(1:50,:),n_state); %descretize according to the rows
a2d=myIntervalDiscretize(data(50:200,:),n_state); %descretize according to the rows
[b,c]=multi_time_series_cat(a1d,a2d);

actualNet=[ 
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
];
t=zeros(1,2);
count=0;
%==================================First iteration: Learning with maxP
maxParentTemp=maxParent;
tic();
[best_net1]=globalMIT_ab_maxP(b,c,alpha,allowSelfLoop,maxParentTemp);
t(1)=toc();

%==================================Finding the nodes that hit maxP
numParents=sum(best_net1);
maxPNodes=find(numParents>=maxParent);
%==================================Second and subsequent iterations:
if(~isempty(maxPNodes))  % if some nodes have hit the max limit
%     fprintf('Nodes hitting maxP limit are: \n');
%     nodeNames(maxPNodes)
  
  while (maxParentTemp<=maxParentLimit)
    %Learning maxPnodes by incrementing maxP
    count=count+1;
    maxParentTemp=maxParentTemp+1;
%     fprintf('Iteration: %d; maxParent: %d\n',count,maxParentTemp);
    [best_net2]=globalMIT_ab_maxP_incr(b,c,alpha,allowSelfLoop,maxParentTemp,maxPNodes);
    t(2)=toc();
    
    numParents=sum(best_net2);
    maxPNodes=find(numParents>=maxParentTemp);    
    
    if(isempty(maxPNodes))  % if some nodes have hit the max limit
        maxParentTemp=maxParentLimit+1;
    end    
    %================================== Combining the networks
    best_net1(:,maxPNodes)=best_net2(:,maxPNodes); %updating the parents for this iteration
  end
end
best_net=best_net1;
fprintf('\nActual Network:\n');
actualNet

fprintf('Learned Network:\n');
best_net

fprintf('Performance Measures:\n True positives, true negatives, false positives, false negatives, precision, recall, f-score, specificity \n');
%tp tn fp fn prec recl fscor spec
M=fnPerformanceMeasure(best_net, actualNet)

fprintf('Time taken in seconds: \n1st iteration: %f\nTotal: %f\n',t(1),t(2));

%creating the network figure using 'createDotGraphic' function
%createDotGraphic(actualNet,nodeNames,'Original SOS network');
%createDotGraphic(best_net,nodeNames,'Learned SOS network');
end
