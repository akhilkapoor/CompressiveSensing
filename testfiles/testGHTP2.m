% This script is used to generate data for the paper !!!!! Change here !!!!!
% J.-L. Bouchot, S. Foucart, P. Hitczenko, "Hard Thresholding Pursuit
% Algorithms: Number of Iterations", 2013

% It analyses the frequency of success and number of iterations of the
% novel algorithm GHTP2. These tests are made for different sparsity levels.


% Author:               Jean-Luc Bouchot, Simon Foucart, Pawel Hitczenko
% Creation date:        07/02/2013
% Modification date:    08/06/2013
% Version:              1
% Copyright:            Math Department, Drexel University, for scholar and
% educational use only


clear all;
close all;
clc;

addpath('./../algorithms/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 200; % size of the measurement vector

NbMat       = 100; % nb of random matrices per sparsities
NbSupport   = 10; % Nb of randomly chosen support per matrix
totTests    = NbMat*NbSupport;


sparsities  = 1:1:floor(2*nbMeasures/3); 

x0          = []; % Where to start the algorithm
MaxNbIter   = nbMeasures;


tolSuccess  = 1e-4; % Tolerance on the relative error for a recovered vector to be consider a success
verbose     = false;

fname = ['GHTP2_gaussVectors_N', num2str(vecSize), '_m', num2str(nbMeasures), '.mat']


%% Generate variables
nbSuccess_ghtp = zeros(length(sparsities),1);

nbIter_ghtp = zeros(length(sparsities),1);

max_nbIter_ghtp = zeros(length(sparsities),1);

nbCorrectIdxPerIter_ghtp = zeros(length(sparsities), vecSize);

min_nbCorrectIdxPerIter_ghtp = sparsities(:) * ones(1,vecSize);

if exist(fname, 'file')
    load(fname);
    firstS = oneS;
    
    nbSuccess_ghtp(firstS) = 0;
    
    nbIter_ghtp(firstS) = 0;
    
    nbCorrectIdxPerIter_ghtp(firstS,:) = 0;
    
    max_nbIter_ghtp(firstS) = 0;
    
    min_nbCorrectIdxPerIter_ghtp(firstS,:) = sparsities(firstS) * ones(1,vecSize);
    
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
            
            
            % Use GHTP
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = ghtp2_gt( y, A, x, { 'initx', x0, 'maxnbiter', MaxNbIter, 'tolres', tolSuccess, 'verbose', verbose } ); 
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_ghtp(thisS) = nbSuccess_ghtp(thisS) + aSuccess;
            if aSuccess
                nbIter_ghtp(thisS) = nbIter_ghtp(thisS) + NbIter;
                max_nbIter_ghtp(thisS) = max(max_nbIter_ghtp(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_ghtp(thisS, oneIter) = nbCorrectIdxPerIter_ghtp(thisS, oneIter) + length(intersect(Ss(oneIter).set,S));
                    min_nbCorrectIdxPerIter_ghtp(thisS, oneIter) = min(min_nbCorrectIdxPerIter_ghtp(thisS, oneIter), length(intersect(Ss(oneIter).set,S)));
                end;
                nbCorrectIdxPerIter_ghtp(thisS,NbIter+1:end) = nbCorrectIdxPerIter_ghtp(thisS,NbIter+1:end) + thisS;
            end;
            
        end;
            
        
    end; % Generate new random matrix
    
    save(fname);
    
end;

%% Average
someSuccessGHTP = find(nbSuccess_ghtp ~= 0);
nbCorrectIdxPerIter_ghtp(someSuccessGHTP,:) = nbCorrectIdxPerIter_ghtp(someSuccessGHTP,:) ./ (nbSuccess_ghtp(someSuccessGHTP,:)*ones(1,size(nbCorrectIdxPerIter_ghtp,2 )));
nbIter_ghtp(someSuccessGHTP) = nbIter_ghtp(someSuccessGHTP) ./ nbSuccess_ghtp(someSuccessGHTP);
nbSuccess_ghtp = 100*nbSuccess_ghtp / totTests;

save(['normalised_', fname]);

