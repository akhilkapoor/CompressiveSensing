%OMP_GT Orthogonal Matching Pursuit with Ground Truth
% OMP_GT(y,A,gt) runs the Orthogonal Matching Pursuit with ground truth
% Usage: [x,S,NormRes,NbIter, Ss, NormRess] = omp_gt(y, A, gt, S0, MaxNbIter, tolRes, verbose)
%   y: measured column vector, 
%   A: measurement matrix, 
%   gt: the original ground truth vector
%   S0: initial index set (optional, default=[]), 
%   MaxNbIter: number of iterations not to be exceeded (optional, default=number of rows of A), 
%   TolRes: tolerance on the relative error for exiting the main loop (optional, default=1e-4)
%   verbose: a boolean true to display some information at each iteration (optional, default=false)
%
%   x: column vector recovered by the algorithm, 
%   S: support of x, 
%   NormRes: the norm of the residual,
%   NbIter: the number of performed iterations
%   Ss: the sequence of active sets selected
%   NormRess: The sequence of norms of residuals
% NOTE: the matrix A is column-normalized before GHTP is run and not
% during this algo
% When using this package, please cite the following paper:
% J.-L. Bouchot, S. Foucart, P. Hitczenko, "Hard Thresholding Pursuit
% Algorithms: Number of Iterations", 2013


% SF (created 25/05/2012, modified 25/05/2012)
% Modifications: JLB, 02/11/13. Added convergence on the sequence of

% Author:               Jean-Luc Bouchot, Simon Foucart, Pawel Hitczenko
% Creation date:        05/25/2013
% Modification date:    06/11/2013
% Version:              1
% Copyright:            Math Department, Drexel University, for scholar and
% educational use only


% This is version developed to use with the ground truth. We get out of the
% loop once all the indices of the original signal have been recovered. 
% the gt argument corresponds to the set of indices in the original signal
function [x,S,NormRes,NbIter, Ss, NormRess] = omp_gt(y, A, gt, S0, MaxNbIter, tolRes, verbose)

[m,N]=size(A);

if nargin < 7 || isempty(verbose)
    verbose = false;
end;

if nargin < 6 || isempty(tolRes)
    tolRes = 1e-4;
end;

if nargin < 5 || isempty(MaxNbIter)
   MaxNbIter=N;
end
if nargin < 4 || isempty(S0)
   S0=[]; 
end

% [A,d]=nzedcol(A);
S=S0;
AS=A(:,S);
xS=AS\y;
res=y-AS*xS;
NormRes=norm(res);
normX = norm(gt);
Sgt = find(gt ~= 0);
NbIter=0;

xS = zeros(N,1);

Ss = [];

while ( NbIter < MaxNbIter && (sum(ismember(Sgt,S)) ~= length(Sgt) ) && (norm(xS-gt) > tolRes*normX) )
    [~,j]=max(abs(A'*res));
    S =sort([S, j]);
    AS=A(:,S);
    xS_=AS\y;
    res=y-AS*xS_;
    xS = zeros(N,1);
    xS(S) = xS_;
    NormRes=norm(res);
    NbIter=NbIter+1;
    NormRess(NbIter) = NormRes;
    Ss(NbIter).set = sort(S);
    
    
    
    if verbose
        disp(['Iter: ', num2str(NbIter), '  ----- NormRes: ', num2str(NormRes), ' ----- Res ratio: ', num2str(NormRes/previousRes), ' -------- Delta: ', num2str(deltaN(NbIter))]);
    end;
    
end

x=zeros(N,1); x(S)=xS_;

end
