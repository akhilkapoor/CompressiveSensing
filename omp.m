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
% Modifications: 
% -> JLB, 02/11/13. Added convergence on the sequence of
% residual, the verbose mode as well as the complete history as output.
% -> JLB, 03/07/2013. Added QR updates in the process to increase the speed
% per iteration.
% 
% NOTE: the matrix A is assumed to be normalized correctly.
function [x,S,NormRes,NbIter, Ss, NormRess, deltaN] = omp(y,A,S0,MaxNbIter,TolRes, verbose, speed)

[m,N]=size(A);

if nargin < 7 || isempty(speed)
    speed = 1; % fast by default
end;
if nargin < 6 || isempty(verbose)
    verbose = false;
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

% % [A,d]=nzedcol(A);
S=S0;
AS=A(:,S);
xS=AS\y;
res=y-AS*xS;
NormRes=norm(res);
NbIter=0;

Q = [];
R = [];
if speed == 1
    [Q,R] = qr(AS);
    linsolveopts.UT = true;
    linsolveopts.LT = false;
end

Ss = [];

previousRes = 2*NormRes;
rn = y;

while ( NbIter < MaxNbIter && NormRes > TolRes )
    [~,j]=max(abs(A'*res));
    S =sort([S, j]);
    AS=A(:,S);
    if speed == 1
        [Q,R] = qrinsert(Q,R, find(S == j),A(:,j));
        xS = linsolve(R, Q'*y, linsolveopts);
    else
        xS=AS\y;
    end
    res=y-AS*xS;
    previousRes = NormRes;
    NormRes=norm(res);
    NbIter=NbIter+1;
    NormRess(NbIter) = NormRes;
    Ss{NbIter} = S;
    
    jidx = find(S == j);
    
    deltaN(NbIter) = xS(jidx)*AS(:,jidx)'*rn;
    rn = res;
    
    
    if verbose
        disp(['Iter: ', num2str(NbIter), '  ----- NormRes: ', num2str(NormRes), ' ----- Res ratio: ', num2str(NormRes/previousRes), ' -------- Delta: ', num2str(deltaN(NbIter))]);
    end;
    
end

% x=zeros(N,1); x(S)=diag(1./d(S))*xS;
x=zeros(N,1); x(S)=xS;

end
