%OMP Orthogonal Matching Pursuit
% OMP(y,A) runs the Orthogonal Matching Pursuit started with x0=0
% Usage: [x,S,NormRes,NbIter] = omp(y,A,S0,MaxNbIter,TolRes)
% y: column vector, A: matrix, S0: initial index set (optional,
% default=[]), MaxNbIter: number of iterations not to be exceeded (optional, default=number of rows of A), TolRes:
% threshold value for the L2-norm of the residual under which the alsogithm
% is stopped (optional, default=1e-4)
% x: column vector, S: support of x, NormRes: the norm of the residual,
% NbIter: the number of performed iterations
% SF (created 25/05/2012, modified 25/05/2012)
% Modifications: JLB, 02/11/13. Added convergence on the sequence of
% residual, the verbose mode as well as the complete history as output.
% NOTE: the matrix A is column-normalized before OMP is run
function [x,S,NormRes,NbIter, Ss, NormRess, deltaN] = omp_sf(y,A,S0,MaxNbIter,TolRes, epsConv, verbose)

[m,N]=size(A);

if nargin < 7 || isempty(verbose)
    verbose = false;
end;
if nargin < 6 || isempty(epsConv)
    epsConv = 1e-2; %this should depend on the noise level in the signal
end;
if nargin < 5 || isempty(TolRes)
   TolRes=1e-4; 
end
if nargin < 4 || isempty(MaxNbIter)
   MaxNbIter=m;
end
if nargin < 3 || isempty(S0)
   S0=[]; 
end

% [A,d]=nzedcol(A);
S=S0;
AS=A(:,S);
xS=AS\y;
res=y-AS*xS;
NormRes=norm(res);
NbIter=0;

Ss = [];

previousRes = 2*NormRes;
rn = y;

while ( NbIter < MaxNbIter && NormRes > TolRes && ~(NormRes/previousRes > 1-epsConv))
    [~,j]=max(abs(A'*res));
    S =sort([S, j]);
    AS=A(:,S);
    xS=AS\y;
    res=y-AS*xS;
    previousRes = NormRes;
    NormRes=norm(res);
    NbIter=NbIter+1;
    NormRess(NbIter) = NormRes;
    Ss{NbIter} = S;
    
    jidx = find(S == j);
%     idx(end)
%     S(idx(end))
%     j
    
    deltaN(NbIter) = xS(jidx)*AS(:,jidx)'*rn;
    rn = res;
    
    
    if verbose
        disp(['Iter: ', num2str(NbIter), '  ----- NormRes: ', num2str(NormRes), ' ----- Res ratio: ', num2str(NormRes/previousRes), ' -------- Delta: ', num2str(deltaN(NbIter))]);
    end;
    
end

x=zeros(N,1); x(S)=xS;

end
