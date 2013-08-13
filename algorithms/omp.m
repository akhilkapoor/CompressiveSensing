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
function [x, S, NormRes, NbIter, Ss, NormRess, deltaN] = omp(y, A, varargin)

[m,N]=size(A);

if size(varargin, 2)
    if iscell(varargin{1})
        [checkArgs, S0, MaxNbIter, TolRes, fast, verbose] = parse_Varargin(m, varargin{1});
    else
        disp('USAGE: ''parameter_name'', parameter_value pairs expected in a cell.');
        disp(['Ignoring arguments: ', varargin]);
        [~, S0, MaxNbIter, TolRes, fast, verbose] = parse_Varargin(m, {});
        checkArgs = false;
        
    end
else
    [checkArgs, S0, MaxNbIter, TolRes, fast, verbose] = parse_Varargin(m, {});
end

if ~checkArgs
%     x = []; S = []; NormRes = 0; NbIter = 0; Ss = []; NormRess = []; deltaN = [];
%     return;
    disp('WARNING:^');
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
if fast == true
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
    if fast == 1
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

function [checkArgs, S0, MaxNbIter, TolRes, Fast, Verbose] = parse_Varargin(m, optional_param)
    
    S0 = [];
    MaxNbIter = m;
    TolRes = 1e-4;
    Verbose = false;
    Fast = false;
    
    [~, s] = size(optional_param);
    if rem(s, 2) ~= 0
        checkArgs = false;
        disp(['Variable argument list is incomplete. Given arguments: ', optional_param{:}]);
%         disp('Given arguments:');
%         disp(optional_param{:});
        disp('USAGE: ''parameter_name'', parameter_value pairs expected in a cell.');
    else
        checkArgs = true;
        for i = 1:2:s
            if strcmpi(optional_param{i}, 'initS')
                S0 = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'MaxNbIter')
                MaxNbIter = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'TolRes')
                TolRes = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'Fast')
                Fast = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'Verbose')
                Verbose = optional_param{i+1};
            else
                checkArgs = false;
                disp(['Unexpected optional parameter: ', optional_param{i}, '. ', mat2str(optional_param{i+1}), ' is being ignored.', optional_param{i}]);
            end
        end
    end
end
