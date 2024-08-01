function [Hz, Hx] = func_split_matrix(H)

    Hz = zeros(size(H));
    Hx = zeros(size(H));

    for i=1:numel(H)
        H_str = num2str(H(i));

        if length(H_str) ~= 2
            Hz(i) = str2num('0');
            Hx(i) = str2num(H_str(2-1));
        else
            Hz(i) = str2num(H_str(1));
            Hx(i) = str2num(H_str(2));
        end
    end
end