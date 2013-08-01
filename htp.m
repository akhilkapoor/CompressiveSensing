%HTP Hard Thresholding Pursuit
% HTP(y,A,s) runs the Hard Thresholding Pursuit started with x0=0
% Usage: [x,S,NormRes,NbIter] = htp(y,A,s,x0,MaxNbIter,TolRes)
% y: column vector, A: matrix, s: integer, x0: initial vector (optional, default=0), 
% MaxNbIter: number of iterations not to be exceeded (optional, default=500), TolRes:
% threshold value for the L2-norm of the residual under which the alsogithm
% is stopped (optional, default=1e-4)
% x: column vector, S: support of x, NormRes: the norm of the residual,
% NbIter: the number of performed iterations
% SF (created 27/05/2012, modified 27/05/2012)
% JLB: Modifications 02/08/2013
% NOTE: the matrix A is column-normalized before HTP is run
function [x,S,NormRes,NbIter, Ss, NormRess] = htp(y,A,s,x0,MaxNbIter,TolRes)

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
    [~,sorted_idx] = sort(abs(u),'descend');
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