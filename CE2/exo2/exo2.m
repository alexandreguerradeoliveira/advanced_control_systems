clc,close all,clear all
load("exo1.mat")
%G2 as the nominal system gives the W2 with the minimum norm, so we should


%% exo2 1
Modulus_margin = 0.5 + 0.01;

W1 = Modulus_margin*tf([1 50],[1 0.000001]);
w1_peak = 20*log10(norm(W1^(-1),inf)) %db peak in high frequencies
W1d = c2d(W1,Ts);
 
% figure 
% bodemag(W1d)

W2 = info_2.W1; % as stated before, take the G2 as nominal system, use matlab's best fit for the model

G = absorbDelay(G2); % absorb delay into the model, matlab does not know how to handle it otherwise
%G = ss(G2); % convert to state space

K = mixsyn(G,W1d,[],W2); % solve the mixed sensitivity problem

W3 = 0.8; % increase W3 => decrease amount of control used
K_nosat = mixsyn(G,W1d,W3,W2); % solve the mixed sensitivity problem

%4
S = feedback(1,G*K);
T = feedback(G*K,1);
U = feedback(K,G);

step(T)
ts = stepinfo(T).SettlingTime % increase the zero in S seems to decrease the ts (at least in LF)

% 5
pzmap(K)
title("pole zero map for controler before reducing order")
K_reduced = minreal(K,0.5);
%T_reduced = feedback(G*K_reduced,1);
figure
pzmap(K_reduced)
title("pole zero map for controler after reducing order")



