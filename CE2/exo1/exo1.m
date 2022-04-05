clc, clear all, close all
%% Params
Ts = 0.04; %[s]
%% CE2 - exo1
%% 1
load data1
Z1 = iddata(y,u,Ts);
Zd1 = detrend(Z1);
G1 = oe(Zd1,[4 6 2]);
G1f = spa(Zd1,200);

%% 2
P = nyquistoptions;
P.ConfidenceRegionDisplaySpacing = 1;
P.Xlim = [-1 1];
P.Ylim = [-.3 1.6];
figure(1);
h = nyquistplot(G1f,P);
axis equal; showConfidence(h,2);

%% 3
figure(2)
h = nyquistplot(G1,P);
axis equal; showConfidence(h,2);

%% 4
%data2
load data2
Z2 = iddata(y,u,Ts);
Zd2 = detrend(Z2);
G2 = oe(Zd2,[4 6 2]);
G1f2 = spa(Zd2,200);

figure(3)
h = nyquistplot(G2,P);
axis equal; showConfidence(h,2);


%data3
load data3
Z3 = iddata(y,u,Ts);
Zd3 = detrend(Z3);
G3 = oe(Zd3,[4 6 2]);
G1f3 = spa(Zd3,200);

figure(4)
h = nyquistplot(G3,P);
axis equal; showConfidence(h,2);

%% 5
Gmm = stack(1,G1,G2,G3);
[Gu_1,info_1] = ucover(Gmm,G1);
W2_g1 = info_1.W1opt;

[Gu_2,info_2] = ucover(Gmm,G2);
W2_g2 = info_2.W1opt;

[Gu_3,info_3] = ucover(Gmm,G3);
W2_g3 = info_3.W1opt;

figure
hold on
bodemag(W2_g1,W2_g2,W2_g3)

% G2 as the nominal system gives the W2 with the minimum norm, so we should
% use it (small w2 => uncertain system is closer to the nominal syst)