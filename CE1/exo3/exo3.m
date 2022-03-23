clc,close all,clear all

k = ureal('k',10,'PlusMinus',[-2 2]);
tau1 = ureal('tau1',2,'PlusMinus',[-0.5 0.5]);
tau2 = ureal('tau2',4,'PlusMinus',[-1 1]);


G1 = tf([k],[tau1 1]);
G2 = tf([1],[tau2 1]);
G = G1*G2;

%nyquist(G);


N1 = 20;
N2 = 200;

% convert to multiplicative uncertanty with 20 samples
g_sample_n1 = usample(G,N1);

order_n1 = 2;

[usys_n1,info_n1] = ucover(g_sample_n1,G.NominalValue,order_n1);

figure
hold on
bodemag(g_sample_n1/G.NominalValue - 1)
bodemag(info_n1.W1)
title("W2 and Ghat/G -1 for 20 samples")
legend("Ghat/G -1 (20 samples)","W2 (second order)")

% convert to multiplicative uncertanty with 200 samples
g_sample_n2 = usample(G,N2);

order_n2 = 2;

[usys_n2,info_n2] = ucover(g_sample_n2,G.NominalValue,order_n2);

figure
hold on
bodemag(g_sample_n2/G.NominalValue - 1)
bodemag(info_n2.W1)
title("W2 and Ghat/G -1 for 200 samples")
legend("Ghat/G -1 (200 samples)","W2 (second order)")




