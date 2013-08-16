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
function [x,S,NormRes,NbIter, Ss, NormRess] = omp_gt(y, A, gt, varargin)

[m,N]=size(A);

if size(varargin, 2)
    if iscell(varargin{1})
        [checkArgs, S0, MaxNbIter, TolRes, verbose] = parse_Varargin(m, varargin{1});
    else
        disp('USAGE: ''parameter_name'', parameter_value pairs expected in a cell.');
        disp(['Ignoring arguments: ', varargin]);
        [~, S0, MaxNbIter, TolRes, verbose] = parse_Varargin(m, {});
        checkArgs = false;
        
    end
else
    [checkArgs, S0, MaxNbIter, TolRes, verbose] = parse_Varargin(m, {});
end

if ~checkArgs
    disp('WARNING:^');
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

while ( NbIter < MaxNbIter && (sum(ismember(Sgt,S)) ~= length(Sgt) ) && (norm(xS-gt) > TolRes*normX) )
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

function [checkArgs, S0, MaxNbIter, TolRes, Verbose] = parse_Varargin(m, optional_param)
    
    S0 = [];
    MaxNbIter = m;
    TolRes = 1e-4;
    Verbose = false;
    
    [~, s] = size(optional_param);
    if rem(s, 2) ~= 0
        checkArgs = false;
        disp('Variable argument list is incomplete.');
        disp('Given arguments:');
        disp(optional_param{:});
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
            elseif strcmpi(optional_param{i}, 'Verbose')
                Verbose = optional_param{i+1};
            else
                checkArgs = false;
                disp(['Unexpected optional parameter: ', optional_param{i}, '. ', mat2str(optional_param{i+1}), ' is being ignored.', optional_param{i}]);
            end
        end
    end
end
