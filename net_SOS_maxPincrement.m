%This file is part of the data and codes used for generating results for the Journal paper:
%*******************************************************************
%Improving gene regulatory network inference using network topology information; 
%A. Nair, M. Chetty, and P. P. Wangikar; Mol. BioSyst., 2015, DOI: 10.1039/C5MB00122F.
%*******************************************************************

%If you find these useful in your work, please cite the above paper.
%(c) 2014-2015 Ajay Nair

%This file is a tutorial on E. coli SOS network reconstruction using maxPincrement algorithm

%Usage:
% net_SOS_maxPincrement()
function []=net_SOS_maxPincrement()
%*********************CONSTANTS REQUIRED
n_state=3;%no of discrete states for the microarray data
alpha=0.999;%significance level for the mutual information test for independance.
allowSelfLoop=0;%allow self regulated link (=1) or not (=0)
maxParent=1; %maximum number of parents a node can have
maxParentLimit=4;
%**************************************
%clc
%**************************************Input data
%the actual E. coli SOS network
actualNet=[
0 0 0 0 0 0 0 0;
1 0 1 0 1 1 1 1;
0 0 0 0 0 0 0 0;
0 1 0 0 0 0 0 0;
0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0];


%network genes
%1-uvrD  2-lexA 3-umuD  4-recA  5-uvrA  6-uvrY  7-ruvA  8-polB
nodeNames=[{'uvrD'},{'lexA'},{'umuD'},{'recA'},{'uvrA'},{'uvrY'},{'ruvA'},{'polB'}];
load data_samples_SOS.mat; %reading the exp data from Ronen-Science2003
%sosPromoterAct 4x8x50 4-experiments; 8-genes; 50-samples
%time  1x50
%*******************************obtaining the experimental data
t1=squeeze(sosPromoterAct(1,:,:));
t2=squeeze(sosPromoterAct(2,:,:));
t3=squeeze(sosPromoterAct(3,:,:));
t4=squeeze(sosPromoterAct(4,:,:));

t1d= myIntervalDiscretize(t1',n_state);%discreatizing
t2d= myIntervalDiscretize(t2',n_state);%discreatizing
t3d= myIntervalDiscretize(t3',n_state);%discreatizing
t4d= myIntervalDiscretize(t4',n_state);%discreatizing
%*******************************
[b,c]=multi_time_series_cat(t1d,t2d,t3d,t4d);



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

