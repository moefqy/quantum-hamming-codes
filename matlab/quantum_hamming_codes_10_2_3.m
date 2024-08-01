% -------------------------------------------------------------------------
% Name         : A. Muh. Mufqi Zuhudi
% NIM          : 1101208451
% Title        : On The Design of Non-CSS Quantum Error Correction Codes 
%                with High Quantum Information
% Email        : moefqy@student.telkomuniversity.ac.id
% Scenario     : Quantum Coding Scheme (Non-CSS) with K=2
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

Hc = [Hc zeros(size(Hc)); zeros(size(Hc)) Hc];

[Hz, Hx] = func_split_matrix(Hc);

% swapped column version
Hc = [01 11 11 01 00 00 00 00 00 00;
      11 11 01 00 00 00 00 00 01 00;
      10 01 01 10 00 00 00 00 00 00;
      01 01 10 00 00 00 00 00 10 00;
      00 00 00 00 01 01 11 11 00 00;
      00 00 00 00 00 11 11 01 00 01;
      00 00 00 00 10 10 01 01 00 00;
      00 00 00 00 00 01 01 10 00 10];
  
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
g5 = func_stabilizer_gen(Hx,Hz,5,'disp');
g6 = func_stabilizer_gen(Hx,Hz,6,'disp');
g7 = func_stabilizer_gen(Hx,Hz,7,'disp');
g8 = func_stabilizer_gen(Hx,Hz,8,'disp');

% define r, N, K
r_H = rank(H);
N = 10;
K = 2;

% perform gauss-jordan elimination to H
H_gauss_jordan = mod(rref(H),2);

% standard form of the encoded Pauli operators X
X_bar1 = func_dynamic_kron(I,Z,Z,I,I,I,I,I,X,I); % input manually
X_bar2 = func_dynamic_kron(I,I,I,I,I,I,Z,Z,I,X); % input manually
 
% standard form of the encoded Pauli operators Z
Z_bar1 = func_dynamic_kron(Z,Z,Z,Z,I,I,I,I,Z,I); % input manually
Z_bar2 = func_dynamic_kron(I,I,I,I,Z,Z,Z,Z,I,Z); % input manually

% transmitter
w1 = eye(length(g1)) + g1;
w2 = eye(length(g2)) + g2;
w3 = eye(length(g3)) + g3;
w4 = eye(length(g4)) + g4;
w5 = eye(length(g5)) + g5;
w6 = eye(length(g6)) + g6;
w7 = eye(length(g7)) + g7;
w8 = eye(length(g8)) + g8;

w_all = w1*w2*w3*w4*w5*w6*w7*w8;

% calculate encoded qubit
ket0000000000 = func_dynamic_kron(ket0,ket0,ket0,ket0,ket0,ket0,ket0,ket0,ket0,ket0);

% 00, 01, 10, 11
ket00L = 1/sqrt(2^r_H) * (w_all * X_bar1^0 * X_bar2^0 * ket0000000000);
ket01L = 1/sqrt(2^r_H) * (w_all * X_bar1^0 * X_bar2^1 * ket0000000000);
ket10L = 1/sqrt(2^r_H) * (w_all * X_bar1^1 * X_bar2^0 * ket0000000000);
ket11L = 1/sqrt(2^r_H) * (w_all * X_bar1^1 * X_bar2^1 * ket0000000000);

% generate transmitted qubit
% b = randi(2,1,1)-1; % random input
b = 11; % input manually

if b == 00
    psiL = ket00L;
elseif b == 01
   psiL = ket01L;
elseif b == 10
   psiL = ket10L;
else 
    psiL = ket11L;
end

% define depolarizing channel
pauli_error = func_dynamic_kron(I,I,I,I,I,I,I,I,I,I);

% apply pauli error to encoded qubit
psiC = pauli_error*psiL;

% measure syndromes
S1 = psiC'*g1*psiC;
S2 = psiC'*g2*psiC;
S3 = psiC'*g3*psiC;
S4 = psiC'*g4*psiC;
S5 = psiC'*g5*psiC;
S6 = psiC'*g6*psiC;
S7 = psiC'*g7*psiC;
S8 = psiC'*g8*psiC;

disp('----------------------------');
func_commute_stabilizer_checker(X_bar1,'disp',g1,g2,g3,g4,g5,g6,g7,g8);
disp('----------------------------');
func_commute_stabilizer_checker(X_bar2,'disp',g1,g2,g3,g4,g5,g6,g7,g8);
disp('----------------------------');
func_commute_stabilizer_checker(Z_bar1,'disp',g1,g2,g3,g4,g5,g6,g7,g8);
disp('----------------------------');
func_commute_stabilizer_checker(Z_bar2,'disp',g1,g2,g3,g4,g5,g6,g7,g8);


S = [S1 S2 S3 S4 S5 S6 S7 S8];
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