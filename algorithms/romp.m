function [x, NbIter] = romp(y, A, varargin)

    % So make the varargin check for only MaxNbIter
    if size(varargin, 2)
        if iscell(varargin{1})
            [checkArgs, MaxNbIter] = parse_Varargin(varargin{1});
        else
            disp('USAGE: ''MaxNbIter'', value expected in a cell.');
            disp(['Ignoring arguments: ', varargin]);
            [~, MaxNbIter] = parse_Varargin({});
            checkArgs = false;
        end
        
    else
        [checkArgs, MaxNbIter] = parse_Varargin({});
    end
    
    if ~checkArgs
        disp('WARNING:^');
    end
    
    [x, NbIter] = romp_Vershynin(MaxNbIter, A, y);
end

function [checkArgs, MaxNbIter] = parse_Varargin(optional_param)

    MaxNbIter = 500;

    [~, s] = size(optional_param);
    if rem(s, 2) ~= 0
        checkArgs = false;
        disp('Variable argument list is incomplete.');
        disp('Given arguments:');
        disp(optional_param{:});
        disp('USAGE: ''MaxNbIter'', value expected.');
    else
        checkArgs = true;
        for i = 1:2:s
            if strcmpi(optional_param{i}, 'MaxNbIter')
                MaxNbIter = optional_param{i+1};
            else
                checkArgs = false;
                disp(['Unexpected optional parameter: ', optional_param{i}, '. ', mat2str(optional_param{i+1}), ' is being ignored.']);
            end
        end
    end
end