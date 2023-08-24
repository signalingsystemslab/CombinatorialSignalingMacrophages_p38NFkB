function [fitI,info1,I,Q]=jacknife_forRoss2014MI(d,c,K,samples,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% JACKNIFE function (adapted from Adelaja et al., 2021) extrapolates mutual information for a set of (multidimensional) inputs,
% X. Actual information calculation is performed by "discrete_continuous_info" (Ross 2014).
%
% d         Labels; d is labels (stimuli/dose etc.) (1 X n_labels)
% c         Metric array;c is your feature matrix that is n_features by n_cells 
% X         Cell array (size=nx1) of sets of individual responses to n different inputs
% K         Kth nearest neighbor (sole parameter passed to KNN algorithm) - not used here
% samples   Number of subsets to measure before computing jacknife extrapolation. Subsets 
%           are randomly drawn to contain between 65%-90% of total samples. The number of 
%           "draws" per sample is computed as 5/(x%)^2.
% varargin  If initialized, forces "verbose mode" (text+graphical output)
%
% fitI      Extrapolated mutual information for entire set
% info1     Mutual information for full set
% I         Mutual information for all subsets
% Q         "Optimal" weights (of all possible inputs) that maximizes mutual information
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<5
    verbose = 0;
else
    verbose = 1;
end

Nsig=length(d);

[info1, Q] = discrete_continuous_info(d,c,3,2); %calculates MI from full set of data 

% #bins, #samples, #draws for jackknife to correct for sample bias        
SampleSize=linspace(0.95,0.6,samples); % set of sample sizes as fractions of total sample size
D= round(5./SampleSize(1:end).^2); % # of draws for each # of samples
Nss=samples; % # of subsets for jackknife

Ik=cell(Nss,1);
for k=1:Nss 
    if verbose; fprintf(['subset' num2str(k)]); end
    info=zeros(D(k),1);
    for ds=1:D(k)    
        if verbose; fprintf('.'); end
            rp=randperm(Nsig);            
            sD = d(:,rp(1:round(SampleSize(k)*Nsig)));
            sC = c(:,rp(1:round(SampleSize(k)*Nsig)));
        info(ds) = discrete_continuous_info(sD,sC,3,2); 
    end
    Ik{k}=info;
end

I=[info1; Ik];

if verbose; fprintf('\n'); end

meanI=zeros(size(I));
stdI=zeros(size(I));
for i=1:length(I)
    meanI(i)=mean(I{i});
    stdI(i)=std(I{i});
end
ss=[1 SampleSize]'; 
pI=polyfit(1./ss,meanI,1);
fitI=pI(2);

% generate plot as jacknife performance control: inverse sample of size size vs mutual information
if verbose
   fig1=figure('position',[50 50 500 500]);
   hold on
   plot(1./ss,meanI,'.')
   plot(1,info1,'.')
   plot([0,2],pI(2)+pI(1)*[0,2],'r')
   xlabel('1/SS')
   ylabel('info (bits)')
end
