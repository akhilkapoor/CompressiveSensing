function [x, Res, NormRes, NormRess, errHist] = CoSaMP(y, A, s, varargin)

    if size(varargin, 2)
        if iscell(varargin{1})
            [checkArgs, errFcn, opts] = parse_Varargin(varargin{1});
        else
            disp('USAGE: ''parameter_name'', parameter_value pairs expected in a cell.');
            disp(['Ignoring arguments: ', varargin]);
            [~, errFcn, opts] = parse_Varargin({});
            checkArgs = false;
        end
        
    else
        [checkArgs, errFcn, opts] = parse_Varargin({});
    end
    
    if ~checkArgs
        disp('WARNING:^');
    end
    
    if ~isfield( opts, 'addK' )
        if s > size(y,1)/3
            error('size of extended support (3*s) cannot be greater than x-dimension of y');
        end
    else
        if s + opts.('addK') > size(y,1)
            error('size of extended support (s+addS) cannot be greater than x-dimension of y');
        end
    end
        
    [x, Res, NormRes, NormRess, errHist] = CoSaMP_Becker(A, y, s, errFcn, opts);

end

function [checkArgs, errFcn, opts] = parse_Varargin(optional_param)

    errFcn = [];
    opts = [];
    MaxNbIter = 500; % maxiter
    % TolRes = 1e-4; % normTol
    % addS = 2*s; % addK
    % TolS = 1e-10; % support_tol
    % HSS = false; % HSS using the HSS Pursuit (see appendix A.2 of Needell/Tropp paper)
    % LSQR_tol % when "A" is a set of function handles, this controls
             % the tolerance in the iterative solver.
    % LSQR_maxNbIter % maximum number of steps in the iterative solver

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
            if strcmpi(optional_param{i}, 'MaxNbIter')
                MaxNbIter = optional_param{i+1};
                opts.('maxiter') = MaxNbIter;
            elseif strcmpi(optional_param{i}, 'tolres')
                TolRes = optional_param{i+1};
                opts.('normTol') = TolRes;
            elseif strcmpi(optional_param{i}, 'adds')
                addS = optional_param{i+1};
                opts.('addK') = addS;
            elseif strcmpi(optional_param{i}, 'tolS')
                TolS = optional_param{i+1};
                opts.('support_tol') = TolS;
            elseif strcmpi(optional_param{i}, 'HSS')
                HSS = optional_param{i+1};
                opts.('HSS') = HSS;
            elseif strcmpi(optional_param{i}, 'LSQR_tol')
                LSQR_tol = optional_param{i+1};
                opts.('LSQR_tol') = LSQR_tol;
            elseif strcmpi(optional_param{i}, 'LSQR_maxNbIter')
                LSQR_maxNbIter = optional_param{i+1};
                opts.('LSQR_maxit') = LSQR_maxNbIter;
            elseif strcmpi(optional_param{i}, 'two_solves')
                two_solves = optional_param{i+1};
                opts.('two_solves') = two_solves;
            elseif strcmpi(optional_param{i}, 'errFcn')
                errFcn = optional_param{i+1};
            else
                checkArgs = false;
                disp('USAGE: ''parameter_name'', parameter_value pairs expected in a cell.');
            end
        end
    end
end