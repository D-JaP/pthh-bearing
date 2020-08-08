function [N1,N2,N3,N4]=hamdang(xi,eta)
    N1=1/4*(1-xi)*(1-eta);
    N2=1/4*(1+xi)*(1-eta);
    N3=1/4*(1+xi)*(1+eta);
    N4=1/4*(1-xi)*(1+eta);
end