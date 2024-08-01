function f = func_dynamic_kron(varargin)

    if nargin < 1
        error('At least one input argument is required.');
    end
    
    f_temp = varargin{1};
    
    for i = 2:nargin
        f_temp = kron(f_temp, varargin{i});
    end

    f = f_temp;
end
