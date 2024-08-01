% -------------------------------------------------------------------------
% Name         : A. Muh. Mufqi Zuhudi
% NIM          : 1101208451
% Title        : On The Design of Non-CSS Quantum Error Correction Codes 
%                with High Quantum Information
% Email        : moefqy@student.telkomuniversity.ac.id
% Scenario     : Quantum Coding Scheme (Non-CSS) with K=1
% -------------------------------------------------------------------------
close all;
clear all;
clc;

% add func script directory
addpath(genpath('func'))

% define qubits
ket0 = [1;0];
ket1 = [0;1];

% define gates
X = [0 1;1 0];
Y = sqrt(-1) * [0 -1;1 0];
Z = [1 0;0 -1];
I = eye(2);

Hdmrd = (1/sqrt(2)) * [1 1;1 -1];
CNOT = [eye(2) zeros(2,2);zeros(2,2) X];

% define H_Hamming
Hc = [01 11 11 01 00;
      11 11 01 00 01;
      10 01 01 10 00;
      01 01 10 00 10];

[Hz, Hx] = func_split_matrix(Hc);

% calculate symplectic inner product (SIP)
SIP = mod((Hx*Hz'+Hz*Hx'),2);

disp('The SIP of Hx and Hz are: ');
disp(SIP);
disp('----------------------------');

% define matrix H using direct sum
H = [Hx Hz];

% define proposed stabilizer
disp('The proposed stabilizers are:')

% define proposed stabilizer X & stabilizer Z
g1 = func_stabilizer_gen(Hx,Hz,1,'disp');
g2 = func_stabilizer_gen(Hx,Hz,2,'disp');
g3 = func_stabilizer_gen(Hx,Hz,3,'disp');
g4 = func_stabilizer_gen(Hx,Hz,4,'disp');

% define r, N, K
r_H = rank(H);
N = 5;
K = 1;

% perform gauss-jordan elimination to H
H_gauss_jordan = mod(rref(H),2);

% standard form of the encoded Pauli operators X
X_bar = func_dynamic_kron(I,Z,Z,I,X); % input manually

% standard form of the encoded Pauli operators Z
Z_bar = func_dynamic_kron(Z,Z,Z,Z,Z);  % input manually

% transmitter
w1 = eye(length(g1)) + g1;
w2 = eye(length(g2)) + g2;
w3 = eye(length(g3)) + g3;
w4 = eye(length(g4)) + g4;
w_all = w1*w2*w3*w4;

% calculate encoded qubit
ket00000 = func_dynamic_kron(ket0,ket0,ket0,ket0,ket0);

ket0L = 1/sqrt(2^r_H) * (w_all*X_bar^0*ket00000);
ket1L = 1/sqrt(2^r_H) * (w_all*X_bar^1*ket00000);

% generate transmitted qubit
% b = randi(2,1,1)-1; % random input
b = 0; % input manually

if b == 0
    psiL = ket0L;
else 
    psiL = ket1L;
end

% define depolarizing channel
pauli_error = func_dynamic_kron(I,I,I,I,I);

% apply pauli error to encoded qubit
psiC = pauli_error*psiL;

% measure syndromes
S1 = psiC'*g1*psiC;
S2 = psiC'*g2*psiC;
S3 = psiC'*g3*psiC;
S4 = psiC'*g4*psiC;

disp('----------------------------');
func_commute_stabilizer_checker(X_bar,'disp',g1,g2,g3,g4);
disp('----------------------------');
func_commute_stabilizer_checker(Z_bar,'disp',g1,g2,g3,g4);

S = [S1 S2 S3 S4];
% normalised_S = (S - min(S)) / (max(S) - min(S));

for i=1:length(S)
    if S(i)>0
    	S(i)=0;
    else
        S(i)=1;
    end
end

disp('----------------------------');
disp('Measured syndrome are:');
disp(S);