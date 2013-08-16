% This script is used to generate data for the paper !!!!! Change here !!!!!
% J.-L. Bouchot, S. Foucart, P. Hitczenko, "Hard Thresholding Pursuit
% Algorithms: Number of Iterations", 2013

% It analyses the frequency of success and number of iterations of the
% non-negative versions of the algorithms. These tests are here for
% linear vectors and for different sparsity levels. 


% Author:               Jean-Luc Bouchot, Simon Foucart, Pawel Hitczenko
% Creation date:        07/06/2013
% Modification date:    08/06/2013
% Version:              1
% Copyright:            Math Department, Drexel University, for scholar and
% educational use only


clear all;
close all;
clc;

addpath('../engine/');

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

fname = ['NN_linearVectors_N', num2str(vecSize), '_m', num2str(nbMeasures), '.mat']


%% Generate variables
nbSuccess_ghtp = zeros(length(sparsities),1);
nbIter_ghtp = zeros(length(sparsities),1);
max_nbIter_ghtp = zeros(length(sparsities),1);
nbCorrectIdxPerIter_ghtp = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_ghtp = sparsities(:) * ones(1,vecSize);

nbSuccess_ghtp2 = zeros(length(sparsities),1);
nbIter_ghtp2 = zeros(length(sparsities),1);
max_nbIter_ghtp2 = zeros(length(sparsities),1);
nbCorrectIdxPerIter_ghtp2 = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_ghtp2 = sparsities(:) * ones(1,vecSize);

nbSuccess_htp = zeros(length(sparsities),1);
nbIter_htp = zeros(length(sparsities),1);
max_nbIter_htp = zeros(length(sparsities),1);
nbCorrectIdxPerIter_htp = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_htp = sparsities(:) * ones(1,vecSize);

nbSuccess_omp = zeros(length(sparsities),1);
nbIter_omp = zeros(length(sparsities),1);
max_nbIter_omp = zeros(length(sparsities),1);
nbCorrectIdxPerIter_omp = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_omp = sparsities(:) * ones(1,vecSize);


% Non-negative algorithms
nbSuccess_ghtpNN = zeros(length(sparsities),1);
nbIter_ghtpNN = zeros(length(sparsities),1);
max_nbIter_ghtpNN = zeros(length(sparsities),1);
nbCorrectIdxPerIter_ghtpNN = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_ghtpNN = sparsities(:) * ones(1,vecSize);

nbSuccess_ghtp2NN = zeros(length(sparsities),1);
nbIter_ghtp2NN = zeros(length(sparsities),1);
max_nbIter_ghtp2NN = zeros(length(sparsities),1);
nbCorrectIdxPerIter_ghtp2NN = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_ghtp2NN = sparsities(:) * ones(1,vecSize);

nbSuccess_htpNN = zeros(length(sparsities),1);
nbIter_htpNN = zeros(length(sparsities),1);
max_nbIter_htpNN = zeros(length(sparsities),1);
nbCorrectIdxPerIter_htpNN = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_htpNN = sparsities(:) * ones(1,vecSize);

nbSuccess_ompNN = zeros(length(sparsities),1);
nbIter_ompNN = zeros(length(sparsities),1);
max_nbIter_ompNN = zeros(length(sparsities),1);
nbCorrectIdxPerIter_ompNN = zeros(length(sparsities), vecSize);
min_nbCorrectIdxPerIter_ompNN = sparsities(:) * ones(1,vecSize);

if exist(fname, 'file')
    %%%%%% HAVE TO ADAPT THIS!!!!!!! %%%%%%%
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
    
	values = 1:-1/thisS:1/thisS;

    
    for oneTrial = 1:NbMat
        
        A = randn( nbMeasures, vecSize )/sqrt(nbMeasures); 
        
        for oneSupp = 1:NbSupport
        
            S = randsample(vecSize, thisS);

            
            x = zeros(vecSize,1);
            x(S) = values;
            
            y = A*x;
            
            
            % Use GHTP2NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = ghtp2_gt( y, A, x, x0, MaxNbIter, tolSuccess, verbose ); 
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_ghtp2(thisS) = nbSuccess_ghtp2(thisS) + aSuccess;
            if aSuccess
                nbIter_ghtp2(thisS) = nbIter_ghtp2(thisS) + NbIter;
                max_nbIter_ghtp2(thisS) = max(max_nbIter_ghtp2(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_ghtp2(thisS, oneIter) = nbCorrectIdxPerIter_ghtp2(thisS, oneIter) + length(intersect(Ss(oneIter).set,S));
                    min_nbCorrectIdxPerIter_ghtp2(thisS, oneIter) = min(min_nbCorrectIdxPerIter_ghtp2(thisS, oneIter), length(intersect(Ss(oneIter).set,S)));
                end;
                nbCorrectIdxPerIter_ghtp2(thisS,NbIter+1:end) = nbCorrectIdxPerIter_ghtp2(thisS,NbIter+1:end) + thisS;
            end;
 
            % Use GHTP_NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = ghtp_gt( y, A, x, x0, MaxNbIter, tolSuccess, verbose );
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
            
            % Use HTP_NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = htp(y,A,thisS,x0,MaxNbIter,tolSuccess); 
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_htp(thisS) = nbSuccess_htp(thisS) + aSuccess;
            if aSuccess
                nbIter_htp(thisS) = nbIter_htp(thisS) + NbIter;
                max_nbIter_htp(thisS) = max(max_nbIter_htp(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_htp(thisS, oneIter) = nbCorrectIdxPerIter_htp(thisS, oneIter) + length(intersect(Ss(oneIter,:),S));
                    min_nbCorrectIdxPerIter_htp(thisS, oneIter) = min(min_nbCorrectIdxPerIter_htp(thisS, oneIter), length(intersect(Ss(oneIter,:),S)));
                end;
                nbCorrectIdxPerIter_htp(thisS,NbIter+1:end) = nbCorrectIdxPerIter_htp(thisS,NbIter+1:end) + thisS;
            end;
            
            % Use OMP_NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = omp_gt(y, A, x, find(x0 ~=0), MaxNbIter, tolSuccess, verbose);
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_omp(thisS) = nbSuccess_omp(thisS) + aSuccess;
            if aSuccess
                nbIter_omp(thisS) = nbIter_omp(thisS) + NbIter;
                max_nbIter_omp(thisS) = max(max_nbIter_omp(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_omp(thisS, oneIter) = nbCorrectIdxPerIter_omp(thisS, oneIter) + length(intersect(Ss(oneIter).set,S));
                    min_nbCorrectIdxPerIter_omp(thisS, oneIter) = min(min_nbCorrectIdxPerIter_omp(thisS, oneIter), length(intersect(Ss(oneIter).set,S)));
                end;
                nbCorrectIdxPerIter_omp(thisS,NbIter+1:end) = nbCorrectIdxPerIter_omp(thisS,NbIter+1:end) + thisS;
            end;
            
            %% Non negative algorithms
            % Use GHTP2NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = ghtp2nn_gt( y, A, x, x0, MaxNbIter, tolSuccess, verbose ); 
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_ghtp2NN(thisS) = nbSuccess_ghtp2NN(thisS) + aSuccess;
            if aSuccess
                nbIter_ghtp2NN(thisS) = nbIter_ghtp2NN(thisS) + NbIter;
                max_nbIter_ghtp2NN(thisS) = max(max_nbIter_ghtp2NN(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_ghtp2NN(thisS, oneIter) = nbCorrectIdxPerIter_ghtp2NN(thisS, oneIter) + length(intersect(Ss(oneIter).set,S));
                    min_nbCorrectIdxPerIter_ghtp2NN(thisS, oneIter) = min(min_nbCorrectIdxPerIter_ghtp2NN(thisS, oneIter), length(intersect(Ss(oneIter).set,S)));
                end;
                nbCorrectIdxPerIter_ghtp2NN(thisS,NbIter+1:end) = nbCorrectIdxPerIter_ghtp2NN(thisS,NbIter+1:end) + thisS;
            end;
 
            % Use GHTP_NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = ghtpnn_gt( y, A, x, x0, MaxNbIter, tolSuccess, verbose );
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_ghtpNN(thisS) = nbSuccess_ghtpNN(thisS) + aSuccess;
            if aSuccess
                nbIter_ghtpNN(thisS) = nbIter_ghtpNN(thisS) + NbIter;
                max_nbIter_ghtpNN(thisS) = max(max_nbIter_ghtpNN(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_ghtpNN(thisS, oneIter) = nbCorrectIdxPerIter_ghtpNN(thisS, oneIter) + length(intersect(Ss(oneIter).set,S));
                    min_nbCorrectIdxPerIter_ghtpNN(thisS, oneIter) = min(min_nbCorrectIdxPerIter_ghtpNN(thisS, oneIter), length(intersect(Ss(oneIter).set,S)));
                end;
                nbCorrectIdxPerIter_ghtpNN(thisS,NbIter+1:end) = nbCorrectIdxPerIter_ghtpNN(thisS,NbIter+1:end) + thisS;
            end;
            
            % Use HTP_NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = htpnn(y,A,thisS,x0,MaxNbIter,tolSuccess); 
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_htpNN(thisS) = nbSuccess_htpNN(thisS) + aSuccess;
            if aSuccess
                nbIter_htpNN(thisS) = nbIter_htpNN(thisS) + NbIter;
                max_nbIter_htpNN(thisS) = max(max_nbIter_htpNN(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_htpNN(thisS, oneIter) = nbCorrectIdxPerIter_htpNN(thisS, oneIter) + length(intersect(Ss(oneIter,:),S));
                    min_nbCorrectIdxPerIter_htpNN(thisS, oneIter) = min(min_nbCorrectIdxPerIter_htpNN(thisS, oneIter), length(intersect(Ss(oneIter,:),S)) );
                end;
                nbCorrectIdxPerIter_htpNN(thisS,NbIter+1:end) = nbCorrectIdxPerIter_htpNN(thisS,NbIter+1:end) + thisS;
            end;
            
            % Use OMP_NN
            [xstar,Sstar,NormRes,NbIter, Ss, NormRess] = ompnn_gt(y, A, x, find(x0 ~= 0), MaxNbIter, tolSuccess, verbose);
            aSuccess = (norm(x-xstar) < tolSuccess*norm(x)); 
            nbSuccess_ompNN(thisS) = nbSuccess_ompNN(thisS) + aSuccess;
            if aSuccess
                nbIter_ompNN(thisS) = nbIter_ompNN(thisS) + NbIter;
                max_nbIter_ompNN(thisS) = max(max_nbIter_ompNN(thisS),NbIter);
                
                for oneIter=1:NbIter
                    nbCorrectIdxPerIter_ompNN(thisS, oneIter) = nbCorrectIdxPerIter_ompNN(thisS, oneIter) + length(intersect(Ss(oneIter).set,S));
                    min_nbCorrectIdxPerIter_ompNN(thisS, oneIter) = min(min_nbCorrectIdxPerIter_ompNN(thisS, oneIter), length(intersect(Ss(oneIter).set,S)));
                end;
                nbCorrectIdxPerIter_ompNN(thisS,NbIter+1:end) = nbCorrectIdxPerIter_ompNN(thisS,NbIter+1:end) + thisS;
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

someSuccessGHTPNN = find(nbSuccess_ghtpNN ~= 0);
nbCorrectIdxPerIter_ghtpNN(someSuccessGHTPNN,:) = nbCorrectIdxPerIter_ghtpNN(someSuccessGHTPNN,:) ./ (nbSuccess_ghtpNN(someSuccessGHTPNN,:)*ones(1,size(nbCorrectIdxPerIter_ghtpNN,2 )));
nbIter_ghtpNN(someSuccessGHTPNN) = nbIter_ghtpNN(someSuccessGHTPNN) ./ nbSuccess_ghtpNN(someSuccessGHTPNN);
nbSuccess_ghtpNN = 100*nbSuccess_ghtpNN / totTests;

someSuccessGHTP2 = find(nbSuccess_ghtp2 ~= 0);
nbCorrectIdxPerIter_ghtp2(someSuccessGHTP2,:) = nbCorrectIdxPerIter_ghtp2(someSuccessGHTP2,:) ./ (nbSuccess_ghtp2(someSuccessGHTP2,:)*ones(1,size(nbCorrectIdxPerIter_ghtp2,2 )));
nbIter_ghtp2(someSuccessGHTP2) = nbIter_ghtp2(someSuccessGHTP2) ./ nbSuccess_ghtp2(someSuccessGHTP2);
nbSuccess_ghtp2 = 100*nbSuccess_ghtp2 / totTests;

someSuccessGHTP2NN = find(nbSuccess_ghtp2NN ~= 0);
nbCorrectIdxPerIter_ghtp2NN(someSuccessGHTP2NN,:) = nbCorrectIdxPerIter_ghtp2NN(someSuccessGHTP2NN,:) ./ (nbSuccess_ghtp2NN(someSuccessGHTP2NN,:)*ones(1,size(nbCorrectIdxPerIter_ghtp2NN,2 )));
nbIter_ghtp2NN(someSuccessGHTP2NN) = nbIter_ghtp2NN(someSuccessGHTP2NN) ./ nbSuccess_ghtp2NN(someSuccessGHTP2NN);
nbSuccess_ghtp2NN = 100*nbSuccess_ghtp2NN / totTests;

someSuccessHTP = find(nbSuccess_htp ~= 0);
nbCorrectIdxPerIter_htp(someSuccessHTP,:) = nbCorrectIdxPerIter_htp(someSuccessHTP,:) ./ (nbSuccess_htp(someSuccessHTP,:)*ones(1,size(nbCorrectIdxPerIter_htp,2 )));
nbIter_htp(someSuccessHTP) = nbIter_htp(someSuccessHTP) ./ nbSuccess_htp(someSuccessHTP);
nbSuccess_htp = 100*nbSuccess_htp / totTests;

someSuccessHTPNN = find(nbSuccess_htpNN ~= 0);
nbCorrectIdxPerIter_htpNN(someSuccessHTPNN,:) = nbCorrectIdxPerIter_htpNN(someSuccessHTPNN,:) ./ (nbSuccess_htpNN(someSuccessHTPNN,:)*ones(1,size(nbCorrectIdxPerIter_htpNN,2 )));
nbIter_htpNN(someSuccessHTPNN) = nbIter_htpNN(someSuccessHTPNN) ./ nbSuccess_htpNN(someSuccessHTPNN);
nbSuccess_htpNN = 100*nbSuccess_htpNN / totTests;

someSuccessOMP = find(nbSuccess_omp ~= 0);
nbCorrectIdxPerIter_omp(someSuccessOMP,:) = nbCorrectIdxPerIter_omp(someSuccessOMP,:) ./ (nbSuccess_omp(someSuccessOMP,:)*ones(1,size(nbCorrectIdxPerIter_omp,2 )));
nbIter_omp(someSuccessOMP) = nbIter_omp(someSuccessOMP) ./ nbSuccess_omp(someSuccessOMP);
nbSuccess_omp = 100*nbSuccess_omp / totTests;

someSuccessOMPNN = find(nbSuccess_ompNN ~= 0);
nbCorrectIdxPerIter_ompNN(someSuccessOMPNN,:) = nbCorrectIdxPerIter_ompNN(someSuccessOMPNN,:) ./ (nbSuccess_ompNN(someSuccessOMPNN,:)*ones(1,size(nbCorrectIdxPerIter_ompNN,2 )));
nbIter_ompNN(someSuccessOMPNN) = nbIter_ompNN(someSuccessOMPNN) ./ nbSuccess_ompNN(someSuccessOMPNN);
nbSuccess_ompNN = 100*nbSuccess_ompNN / totTests;


save(['normalised_', fname]);

%% Finish normalization

%% Some first plots
figure; hold on; plot(sparsities, nbSuccess_ghtp, 'r'); plot(sparsities, nbSuccess_ghtpNN, 'b'); legend('GHTP','GHTPNN');
figure; hold on; plot(sparsities, nbSuccess_ghtp2, 'r'); plot(sparsities, nbSuccess_ghtp2NN, 'b'); legend('GHTP2','GHTP2NN');
figure; hold on; plot(sparsities, nbSuccess_htp, 'r'); plot(sparsities, nbSuccess_htpNN, 'b'); legend('HTP','HTPNN');
figure; hold on; plot(sparsities, nbSuccess_omp, 'r'); plot(sparsities, nbSuccess_ompNN, 'b'); legend('OMP','OMPNN');

