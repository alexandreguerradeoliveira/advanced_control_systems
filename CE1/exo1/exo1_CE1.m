clc,close all,clear all
%% exo 1.1
%1.1.1
a = 1;
b = 0.1;
c = 16;

delta = (b^2)-4*a*c;

s1 = (-b+sqrt(delta))/(2*a);
s2 = (-b-sqrt(delta))/(2*a);
s3 = (b+sqrt(delta))/(2*a);
s4 = (b-sqrt(delta))/(2*a);

Res_s1 = (10-2*s1)*(10+2*s1)/((s1-s2)*(s1-s3)*(s1-s4));
Res_s2 = (10-2*s2)*(10+2*s2)/((s2-s1)*(s2-s3)*(s2-s4));
Res_s3 = (10-2*s3)*(10+2*s3)/((s3-s1)*(s3-s2)*(s3-s4));
Res_s4 = (10-2*s4)*(10+2*s4)/((s4-s1)*(s4-s2)*(s4-s3));

norm2_1 = sqrt(Res_s1+Res_s2)

%1.1.2
fun = @(om) (1/(2*pi))*(abs((10-2*(1i*om))./( (((1i*om)).^2)+0.1*(1i*om)+16  )).^2);

norm_2_2 = sqrt(integral(fun,-Inf,Inf))

%1.1.3

G = tf([-2 10],[1 0.1 16]);
[impulse_g,impulse_tt] = impulse(G);

norm_2_3 = sqrt(trapz(impulse_tt,impulse_g.*impulse_g))

%1.1.4
A = [0,1;-16,-0.1];
B = [0;1];
C = [10,-2];

L = are(A',zeros(2,2),B*(B'));

norm_2_4 = sqrt(trace(C*L*(C')))

%1.1.5

norm_2_5 = norm(G,2)

%% exo 1.2
%1.2.1
fun2 = @(om) abs((10-2*(1i*om))./( (((1i*om)).^2)+0.1*(1i*om)+16  ));

figure
bode(G) % we can see that omega_N = 100 [rad/s] is fine because the gain is already very small.

omega_grid = linspace(0,100,10000);
freq_response = fun2(omega_grid);

norm_inf_1 = max(freq_response)

%1.2.2
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

norm_inf_2 = gam
%1.2.3
norm_inf_3 = norm(G,Inf)
