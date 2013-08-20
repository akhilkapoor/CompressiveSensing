function [x, NbIter] = romp(y, A, s, varargin)

    if size(varargin, 2)
        if iscell(varargin{1})
            [checkArgs, MaxNbIter, TolRes] = parse_Varargin(varargin{1});
        else
            disp('USAGE: ''parameter_name'', paramenter_value pairs expected in a cell.');
            disp(['Ignoring arguments: ', varargin]);
            [~, MaxNbIter, TolRes] = parse_Varargin({});
            checkArgs = false;
        end
    else
        [checkArgs, MaxNbIter, TolRes] = parse_Varargin({});
    end
    
    if ~checkArgs
        disp('WARNING:^');
    end
    
    [x, NbIter] = romp_Vershynin(s, A, y, MaxNbIter, TolRes);
end

function [checkArgs, MaxNbIter, TolRes] = parse_Varargin(optional_param)

    MaxNbIter = 500;
    TolRes = 1e-4;

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
            if strcmpi(optional_param{i}, 'MaxNbIter')
                MaxNbIter = optional_param{i+1};
            elseif strcmpi(optional_param{i}, 'TolRes')
                TolRes = optional_param{i+1};
            else
                checkArgs = false;
                disp(['Unexpected optional parameter: ', optional_param{i}, '. ', mat2str(optional_param{i+1}), ' is being ignored.']);
            end
        end
    end
end