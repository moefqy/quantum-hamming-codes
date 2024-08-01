function f = func_stabilizer_gen(Hx,Hz,gn, mode)
    
    % define gates
    X = [0 1;1 0];
    Y = sqrt(-1) * [0 -1;1 0];
    Z = [1 0;0 -1];
    I = eye(2);

    % define matrix H 
    H = [Hx Hz];
    
    g_proposed = '';
    g = [];
    
    for i=1:size(H,2)/2
        if i==1
            if Hx(gn,i) == 1 && Hz(gn,i) == 1 
                g_proposed = strcat(g_proposed, 'Y');
                g = Y;
            elseif Hx(gn,i) == 1
                g_proposed = strcat(g_proposed, 'X');
                g = X;
            elseif Hz(gn,i) == 1
                g_proposed = strcat(g_proposed, 'Z');
                g = Z;
            else
                g_proposed = strcat(g_proposed, 'I');
                g = I;
            end
        else
            if Hx(gn,i) == 1 && Hz(gn,i) == 1 
                g_proposed = strcat(g_proposed, 'Y');
                g = kron(g,Y);
            elseif Hx(gn,i) == 1
                g_proposed = strcat(g_proposed, 'X');
                g = kron(g,X);
            elseif Hz(gn,i) == 1
                g_proposed = strcat(g_proposed, 'Z');
                g = kron(g,Z);
            else
                g_proposed = strcat(g_proposed, 'I');
                g = kron(g,I);
            end
        end
    end
    
    switch mode
        case 'disp'
            disp(['g', num2str(gn), ' = ', g_proposed])
            f = g;
        case 'nodisp'
            f = g;
end