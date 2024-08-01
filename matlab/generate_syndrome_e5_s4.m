% -------------------------------------------------------------------------
% Name         : A. Muh. Mufqi Zuhudi
% NIM          : 1101208451
% Title        : On The Design of Non-CSS Quantum Error Correction Codes 
%                with High Quantum Information
% Email        : moefqy@student.telkomuniversity.ac.id
% Scenario     : Measure syndrome of possible error for k=1
% -------------------------------------------------------------------------

% input pauli list manually
pauli_list = ['IIIII';
              'XIIII';
              'IXIII';
              'IIXII';
              'IIIXI';
              'IIIIX';
              'YIIII';
              'IYIII';
              'IIYII';
              'IIIYI';
              'IIIIY';
              'ZIIII';
              'IZIII';
              'IIZII';
              'IIIZI';
              'IIIIZ'];

S_all = [];

for ii=1:size(pauli_list,1)
    
    pauli_str = pauli_list(ii,:);
    
    % convert each character in the pauli string to the corresponding matrix operator
    pauli_ops = {};
    for jj = 1:length(pauli_str)
        switch pauli_str(jj)
            case 'I'
                pauli_ops{end+1} = eye(2);
            case 'X'
                pauli_ops{end+1} = [0 1;1 0];
            case 'Y'
                pauli_ops{end+1} = [0 -1i;1i 0];
            case 'Z'
                pauli_ops{end+1} = [1 0;0 -1];
        end
    end
    
    % calculate the Kronecker product of all the matrix operators
    pauli_error = func_dynamic_kron(pauli_ops{:});
    
    psiC = pauli_error*psiL;
    
    % measure syndromes
    S1 = psiC'*g1*psiC;
    S2 = psiC'*g2*psiC;
    S3 = psiC'*g3*psiC;
    S4 = psiC'*g4*psiC;
    
    S = [S1 S2 S3 S4];

    for kk=1:length(S)
        if S(kk)>0
            S(kk)=0;
        else
            S(kk)=1;
        end
    end
    
    S_all = cat(1,S_all, S);

    disp(['Measured syndrome of ', pauli_str, ' are: ', num2str(S)]);
end

disp('----------------------------');
disp(['Unique syndrome are ', num2str(size(unique(S_all,'rows'),1)), ' out of ', num2str(size(pauli_list,1))])
