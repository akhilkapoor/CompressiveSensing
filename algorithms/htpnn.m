%HTPNN Hard Thresholding Pursuit - Non-Negative
% HTPNN(y,A,s) runs the Hard Thresholding Pursuit started with x0=0
% Usage: [x,S,NormRes,NbIter, Ss, NormRess] = htpnn(y,A,s,x0,MaxNbIter,TolRes)
%   y: column vector: measurements, 
%   A: measurement/sensing matrix, 
%   s: sparsity level
%   x0: initial vector (optional, default=0), 
%   MaxNbIter: number of iterations not to be exceeded (optional, default=500), 
%   tolRes: tolerance on the relative error (optional, default 1e-4)
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
function [x,S,NormRes,NbIter, Ss, NormRess] = htpnn(y,A,s,x0,MaxNbIter,TolRes)

[~,N]=size(A);
if nargin < 6 || isempty(TolRes)
   TolRes=1e-4; 
end
if nargin < 5 || isempty(MaxNbIter)
   MaxNbIter=500;
end
if nargin < 4 || isempty(x0)
   x0=zeros(N,1); 
end

x = x0;
res = y-A*x;
NormRes = norm(res);
NbIter=0;

[~, idx] = sort(abs(x0));
S = idx(1:s);
Snew = idx(end-s+1:end); % Assign something completely different to make sure wwe get at least once in the for loop

Ss = [];

while ( (sum(S==Snew)<s) && NbIter < MaxNbIter && NormRes > TolRes )
    u = x+A'*res;
    [~,sorted_idx] = sort(u,'descend');
    S = sorted_idx(1:s);
    x=zeros(N,1); 
    x(S)=A(:,S)\y; % Here some optimization can be done based on QR decomposition and updates
    res = y-A*x;
    NormRes=norm(res);
    NbIter = NbIter+1;
    Ss = [Ss; sort(S(:))'];
    NormRess(NbIter) = NormRes;
    aux = Snew;
    Snew = S;
    S = aux;
%     disp(['NbIter in HTP:', num2str(NbIter)])
%     sum(S==Snew)
%     pause
end

end