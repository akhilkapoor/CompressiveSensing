% This script is a modified version of testGaussianVectors.m to work only
% with "Hard Thresholding Pursuit" and "Orthogonal Matching Pursuit"

% It considers frequency of success, number of iterations, number of
% correct indices, norm of the output, etc... for the case of gaussian
% vectors. OMP, HTP and GHTP are compared. 
%
% Two .mat files will be created, one containing raw datas and the second
% one containing normalized data that will be used to graph (see graph*.m
% files) 
% NOTE: The files created contain the size of the measurement matrix in
% their names. Changing the size will create a new file, but changing any
% other parameter will reuse the exact same size. Make sure to rename the
% matlab files if needed.

% Author:               Jean-Luc Bouchot, Simon Foucart, Pawel Hitczenko
% Creation date:        04/29/2013
% Modification date:    06/11/2013
% Version:              1
% Copyright:            Math Department, Drexel University, for scholar and
% educational use only


clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 200; % size of the measurement vector

NbMat       = 100; % nb of random matrices per sparsities
NbSupport   = 10; % Nb of randomly chosen support per matrix
totTests    = NbMat*NbSupport;


sparsities  = 1:1:90; 

x0          = []; % Where to start the algorithm
MaxNbIter   = nbMeasures;


tolSuccess  = 1e-4; % Tolerance on the relative error for a recovered vector to be consider a success
verbose     = false;

fname = ['gaussVectors_N', num2str(vecSize), '_m', num2str(nbMeasures), '.mat']


%% Generate variables
nbSuccess_omp = zeros(length(sparsities),1);
nbSuccess_htp = zeros(length(sparsities),1);

nbIter_htp = zeros(length(sparsities),1);
nbIter_omp = zeros(length(sparsities),1);

max_nbIter_htp = zeros(length(sparsities),1);
max_nbIter_omp = zeros(length(sparsities),1);

nbCorrectIdxPerIter_omp = zeros(length(sparsities), vecSize);
nbCorrectIdxPerIter_htp = zeros(length(sparsities), vecSize);

min_nbCorrectIdxPerIter_omp = sparsities(:) * ones(1,vecSize);
min_nbCorrectIdxPerIter_htp = sparsities(:) * ones(1,vecSize);

if exist(fname, 'file')
    load(fname);
    firstS = oneS;
    
    nbSuccess_omp(firstS) = 0;
    nbSuccess_htp(firstS) = 0;
    
    nbIter_htp(firstS) = 0;
    nbIter_omp(firstS) = 0;
    
    nbCorrectIdxPerIter_omp(firstS,:) = 0;
    nbCorrectIdxPerIter_htp(firstS,:) = 0;
    
    max_nbIter_htp(firstS) = 0;
    max_nbIter_omp(firstS) = 0;
    
    min_nbCorrectIdxPerIter_omp(firstS,:) = sparsities(firstS) * ones(1,vecSize);
    min_nbCorrectIdxPerIter_htp(firstS,:) = sparsities(firstS) * ones(1,vecSize);
    
else
    firstS = 1;
end;

%% Launch tests
for oneS=firstS:length(sparsities)
    thisS = sparsities(oneS)
    S = 1:thisS;
    
    for oneTrial = 1:NbMat
        
        A = randn( nbMeasures, vecSize )/sqrt(nbMeasures); 
        
        for oneSupp = 1:NbSupport
            
            x = zeros(vecSize,1);
            x(S) = randn(thisS,1);
            
            y = A*x;
            
            % Use HTP
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = htp(y,A,thisS,x0,MaxNbIter,0);  
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_htp(thisS) = nbSuccess_htp(thisS) + aSuccess;
            
            if aSuccess
                nbIter_htp(thisS) = nbIter_htp(thisS) + NbIter;
                max_nbIter_htp(thisS) = max(max_nbIter_htp(thisS), NbIter);
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_htp(thisS, oneIter) = nbCorrectIdxPerIter_htp(thisS, oneIter) + length(intersect(Ss(oneIter,:),S));
                    min_nbCorrectIdxPerIter_htp(thisS, oneIter) = min(min_nbCorrectIdxPerIter_htp(thisS, oneIter), length(intersect(Ss(oneIter,:),S)));
                end;
                nbCorrectIdxPerIter_htp(thisS,NbIter+1:end) = nbCorrectIdxPerIter_htp(thisS,NbIter+1:end) + thisS;
            end;
            
            % Use OMP
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = omp(y,A, x, find(x0 ~=0),MaxNbIter, tolSuccess, verbose); 
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_omp(thisS) = nbSuccess_omp(thisS) + aSuccess;
            
            if aSuccess
                nbIter_omp(thisS) = nbIter_omp(thisS) + NbIter;
                max_nbIter_omp(thisS) = max(max_nbIter_omp(thisS), NbIter);
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_omp(thisS, oneIter) = nbCorrectIdxPerIter_omp(thisS, oneIter) + length(intersect(Ss(oneIter).set,S));
                    min_nbCorrectIdxPerIter_omp(thisS, oneIter) = min(min_nbCorrectIdxPerIter_omp(thisS, oneIter), length(intersect(Ss(oneIter).set,S)));
                end;
                nbCorrectIdxPerIter_omp(thisS,NbIter+1:end) = nbCorrectIdxPerIter_omp(thisS,NbIter+1:end) + thisS;
            end;
            
        end;
            
        
    end; % Generate new random matrix
    
    save(fname);
    
end;

%% Average
someSuccessHTP = find(nbSuccess_htp ~= 0);
someSuccessOMP = find(nbSuccess_omp ~= 0);

nbCorrectIdxPerIter_htp(someSuccessHTP,:) = nbCorrectIdxPerIter_htp(someSuccessHTP,:) ./ (nbSuccess_htp(someSuccessHTP,:)*ones(1,size(nbCorrectIdxPerIter_htp,2 )));
nbCorrectIdxPerIter_omp(someSuccessOMP,:) = nbCorrectIdxPerIter_omp(someSuccessOMP,:) ./ (nbSuccess_omp(someSuccessOMP,:)*ones(1,size(nbCorrectIdxPerIter_omp,2 )));

nbIter_htp(someSuccessHTP) = nbIter_htp(someSuccessHTP) ./ nbSuccess_htp(someSuccessHTP);
nbIter_omp(someSuccessOMP) = nbIter_omp(someSuccessOMP) ./ nbSuccess_omp(someSuccessOMP);

nbSuccess_htp = 100*nbSuccess_htp / totTests;
nbSuccess_omp = 100*nbSuccess_omp / totTests;

save(['normalised_', fname]);
