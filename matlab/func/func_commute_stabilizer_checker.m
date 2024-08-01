function f = func_commute_stabilizer_checker(XZ_bar, mode, varargin)

    bar_name = inputname(1);
    f = 0;
    switch mode
        case 'disp'
            for i=1:size(varargin,2)
                if varargin{i}*XZ_bar == XZ_bar*varargin{i}
                    disp([bar_name, ' and g', num2str(i), ' is commute']) % comment to display
                    f = f+1;
                else
                    disp([bar_name, ' and g', num2str(i), ' is NOT commute']) % comment to display
                    f = f;
                end
            end
        case 'nodisp'
            for i=1:size(varargin,2)
                if varargin{i}*XZ_bar == XZ_bar*varargin{i}
                    f = f+1;
                else
                    f = f;
                end
            end
    end
end
