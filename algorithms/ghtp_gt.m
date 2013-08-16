%GHTP_GT Graded Hard Thresholding Pursuit (with ground truth)
% GHTP(y,A,gt) runs the Hard Thresholding Pursuit started with x0=0 and with
% an increasing size of the support
% Usage: [ x,S,NormRes,NbIter, Ss, NormRess ] = ghtp_gt( y, A, gt, x0, MaxNbIter, tolRes, verboseMode )
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
% When using this package, please cite the following paper:
% J.-L. Bouchot, S. Foucart, P. Hitczenko, "Hard Thresholding Pursuit
% Algorithms: Number of Iterations", 2013

% Author:               Jean-Luc Bouchot, Simon Foucart, Pawel Hitczenko
% Creation date:        05/25/2013
% Modification date:    06/11/2013
% Version:              1
% Copyright:            Math Department, Drexel University, for scholar and
% educational use only

function [ x,S,NormRes,NbIter, Ss, NormRess ] = ghtp_gt( y, A, gt, varargin )

[~,N]=size(A);

if size(varargin, 2)
    [checkArgs, x0, MaxNbIter, TolRes, verbose] = parse_Varargin(N, varargin{:});
else
    [checkArgs, x0, MaxNbIter, TolRes, verbose] = parse_Varargin(N, {});
end

if ~checkArgs
    disp('WARNING:^');
end

x = x0;
S = find(x ~= 0);
res = y-A*x;
NormRes = norm(res);
NbIter=0;

normX = norm(gt);
Sgt = find(gt ~= 0);

while ( NbIter < MaxNbIter && (sum(ismember(Sgt,S)) ~= length(Sgt) ) && (norm(x-gt) > TolRes*normX) )
    u = x+A'*res;
    [~,sorted_idx] = sort(abs(u),'descend');
    S = sorted_idx(1:NbIter+1);
    x=zeros(N,1); 
    x(S)=A(:,S)\y; % Improvement possible by updating the QR decomposition
    res = y-A*x;
    previousRes = NormRes;
    NormRes=norm(res);
    NbIter = NbIter+1;
    if verbose
        disp(['Iter: ', num2str(NbIter), '  ----- NormRes: ', num2str(NormRes), ' ----- Res ratio: ', num2str(NormRes/previousRes)]);
    end;
    Ss(NbIter).set = sort(S);
    NormRess(NbIter) = NormRes;
end

end

function [checkArgs, x0, MaxNbIter, TolRes, Verbose] = parse_Varargin(N, optional_param)
    
    x0=zeros(N,1);
    MaxNbIter = 500;
    TolRes = 1e-4;
    Verbose = false;
    
    [~, s] = size(optional_param);
    if rem(s, 2) ~= 0
        checkArgs = false;
        disp('Variable argument list is incomplete.');
        disp('Given arguments:');
        disp(optional_param{:});
        disp('USAGE: ''parameter_name'', parameter_value pairs expected.');
    else
        checkArgs = true;
        for i = 1:2:s
            if strcmpi(optional_param{i}, 'initX')
                x0 = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'MaxNbIter')
                MaxNbIter = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'TolRes')
                TolRes = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'Verbose')
                Verbose = optional_param{i+1};
            else
                checkArgs = false;
                disp(['Unexpected optional parameter: ', optional_param{i}, '. ', mat2str(optional_param{i+1}), ' is being ignored.']);
            end
        end
    end
end
