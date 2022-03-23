clc,close all,clear all
data = importdata('G_mimo.mat');


%% continous time state space model
A = data.A;
B = data.B;
C = data.C;
D = data.D;

sys_cts = ss(A,B,C,D);


%% exercice 1.2.1

% 1. 2-norm with frequency responce
sigma(sys_cts) % plot the sigma to know where to stop and end the grid in omega

N = 50000;
omega_grid = logspace(-4,4,N);
G_grided = freqresp(sys_cts,omega_grid);

%integrate_function = pagemtimes(pagetranspose(G_grided),G_grided);

%integrate_function_trace = zeros(1,N);
for k = 1:N
    integrate_function(:,:,k)= (G_grided(:,:,k)')*(G_grided(:,:,k));
    integrate_function_trace(k) = trace(integrate_function(:,:,k));
end
 % function is symetric so integrate from 0 to inf then multiply by 2
estimated_integral = 2*trapz(omega_grid,integrate_function_trace);

norm_2_1 = sqrt((1/(2*pi))*(estimated_integral))

% 2. 2-norm with state space method

L = are(A',zeros(size(A)),B*(B'));

norm_2_q2 = sqrt(trace(C*L*(C')))


% 3. 2-norm with matlab


norm_2_q3 = norm(sys_cts,2)

%% exercice 1.2.2

% 1. inf-norm with frequency response

[sig_frequency_response,~] = sigma(sys_cts);

max_sig_frequency_response = max(sig_frequency_response,[],1);

norm_inf_q1 = max(max_sig_frequency_response)

% 2. inf-norm with the bonded real lemma
H = @(gam) [A,(gam^(-2))*B*(B');-C'*C,-A'];

tol_gam = 10^(-10);
tol_eigs = 10^(-10);

gam1 = 0.1;
gam2 = 100;

while abs((gam1-gam2)/gam1)>tol_gam
    % solve for the eigen values
    
    gam = (gam2+gam1)/2;
    eig_gam = eig(H(gam));
    
    up = 0; % up = 0 <=> H has no eigenvalue in the imag axis for the current gamma
    
    for k = 1:length(eig_gam)
       if (abs(real(eig_gam(k)))<=tol_eigs) 
          up = 1; % H has a eigenvalue on the imaginary axis (or very close) for this gamma
       end
    end
    
    if (up==0)  
        gam2 = gam; % if H has no an eigenvalue on the imaginary axis set gam2 = gam
    else
        gam1 = gam; % if H has an eigenvalue on the imaginary axis set gam1 = gam
    end
end

norm_inf_q2 = gam

% 3. inf-norm with matlab command

norm_inf_q3 = norm(sys_cts,inf)