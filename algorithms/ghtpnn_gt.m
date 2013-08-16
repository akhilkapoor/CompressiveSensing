%GHTPNN_GT Graded Hard Thresholding Pursuit Non-Negative (with ground truth)
% GHTPNN_GT(y,A,gt) runs the Hard Thresholding Pursuit started with x0=0 and with
% an increasing size of the support
% Usage: [ x,S,NormRes,NbIter, Ss, NormRess ] = ghtpnn_gt( y, A, gt, x0, MaxNbIter, tolRes, verboseMode )
%   y: column vector, 
%   A: matrix, 
%   gt: ground truth for the vector to be recovered, 
%   x0: initial vector (optional, default=0), 
%   MaxNbIter: number of iterations not to be exceeded (optional, default=500), 
%   tolRes: tolerance on the relative error (optional, default 1e-4)
%   verboseMode: A boolean controling the amount of information displayed (optional, default=false)
% 
%   x: column vector recovered by the algorithm
%   S: support of x, 
%   NormRes: the norm of the residual when exiting the main loop
%   NbIter: the number of performed iterations
%   Ss: the sequence of active sets selected
%   NormRess: The sequence of norms of residuals
% NOTE: the matrix A is column-normalized before GHTP is run and not
% during this algo
% When using this package, please cite the following paper: !!!!! Change here !!!!!
% J.-L. Bouchot, S. Foucart, P. Hitczenko, "Hard Thresholding Pursuit
% Algorithms: Number of Iterations", 2013

% Author:               Jean-Luc Bouchot, Simon Foucart, Pawel Hitczenko
% Creation date:        07/02/2013
% Modification date:    07/02/2013
% Version:              1
% Copyright:            Math Department, Drexel University, for scholar and
% educational use only

function [ x,S,NormRes,NbIter, Ss, NormRess ] = ghtpnn_gt( y, A, gt, x0, MaxNbIter, tolRes, verboseMode )

[~,N]=size(A);

if nargin < 7 || isempty(verboseMode)
    verboseMode = false;
end;

if nargin < 6 || isempty(tolRes)
	tolRes = 1e-4;
end;


if nargin < 5 || isempty(MaxNbIter)
   MaxNbIter=500;
end

if nargin < 4 || isempty(x0)
   x0=zeros(N,1); 
end

x = x0;
S = find(x ~= 0);
res = y-A*x;
NormRes = norm(res);
NbIter=0;

normX = norm(gt);
Sgt = find(gt ~= 0);

while ( NbIter < MaxNbIter && (sum(ismember(Sgt,S)) ~= length(Sgt) ) && (norm(x-gt) > tolRes*normX) )
    u = x+A'*res;
    [~,sorted_idx] = sort(u,'descend');
    S = sorted_idx(1:NbIter+1);
    x=zeros(N,1); 
    x(S)=A(:,S)\y; % Improvement possible by updating the QR decomposition
    res = y-A*x;
    previousRes = NormRes;
    NormRes=norm(res);
    NbIter = NbIter+1;
    if verboseMode
        disp(['Iter: ', num2str(NbIter), '  ----- NormRes: ', num2str(NormRes), ' ----- Res ratio: ', num2str(NormRes/previousRes)]);
    end;
    Ss(NbIter).set = sort(S);
    NormRess(NbIter) = NormRes;
end

end

