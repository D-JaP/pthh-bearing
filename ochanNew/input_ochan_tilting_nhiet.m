global beta r1 r2 nuy omega theta_p r_p X_p Y_p kchia_m kchia_n kchia_q T0 Co_h Co_k t_pad nuy0
% beta=20/360*2*pi;
% r1=1.105;
% r2=1.690;
% nuy=0.004043;%do n0hot dong luc hoc (N/m^2.s)
% omega=107/60*(2*pi);%toc do goc (rad/s)
% theta_p=11.415/180*pi
% r_p=1.4175
% X_p=r_p*sin(theta_p);
% Y_p=r_p*cos(theta_p);

% beta=55/360*2*pi;
% r1=0.02;
% r2=0.05;
% nuy=0.039;%do nhot dong luc hoc (N/m^2.s)
% omega=3000/60*(2*pi);%toc do goc (rad/s)

% beta=50/360*2*pi;
% r1=0.05715;
% r2=0.1143;
% nuy=0.042646;%do nhot dong luc hoc (N/m^2.s)
% omega=3000/60*(2*pi);%toc do goc (rad/s)
% theta_p=30/180*pi
% r_p=0.085725
% X_p=r_p*sin(theta_p);
% Y_p=r_p*cos(theta_p);

kchia_m =30;    kchia_n=30; kchia_q=10;
beta=28/360*2*pi;
r1=0.1875;
r2=0.3225;
nuy0=0.0252;%do nhot dong luc hoc (N/m^2.s)
omega=3000/60*(2*pi);%toc do goc (rad/s)
theta_p=17.38/180*pi;
r_p=0.255;
X_p=r_p*sin(theta_p);
Y_p=r_p*cos(theta_p);


% beta=30/360*2*pi;
% r1=0.0285;
% r2=0.045;
% nuy0=0.0220;%do nhot dong luc hoc (N/m^2.s)
% omega=500/60*(2*pi);%toc do goc (rad/s)
% theta_p=1/180*pi;
% r_p=0.0367;
% X_p=r_p*sin(theta_p);
% Y_p=r_p*cos(theta_p);


k_conduct=0.13; 
delta=855;%khoi luong
Cp=2090; %capacity



T0=30;
t_pad=0.1;
Co_k=1;
Co_h=1;
T_alpha=1;
%khong thu nguyen
% nuy0=1;%%%%%%%%%%%%%%%%%%%%%%%%inlet oil viscosity
