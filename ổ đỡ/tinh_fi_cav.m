function [f,gamma_cav]=tinh_fi_cav(epsilon)
gamma_cav=4;
k=1;deltaf=1;
while (k<=30&&abs(deltaf)>=1e-5)
    f=epsilon*(gamma_cav-sin(gamma_cav)*cos(gamma_cav))-2*sin(gamma_cav)+2*gamma_cav*cos(gamma_cav);
    df=epsilon*(1-cos(gamma_cav)^2+sin(gamma_cav)^2)-2*gamma_cav*sin(gamma_cav);
    deltaf=f/df;
    gamma_cav=gamma_cav-deltaf;
    k=k+1;
end

fi_cav=2*pi-acos((epsilon-cos(gamma_cav))/(epsilon*cos(gamma_cav)-1));
f=fi_cav;