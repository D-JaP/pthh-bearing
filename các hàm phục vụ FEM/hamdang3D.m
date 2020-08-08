function [N1,N2,N3,N4,N5,N6,N7,N8]=hamdang3D(r,s,t)
    N1=1/8*(1-r)*(1-s)*(1-t);
    N2=1/8*(1+r)*(1-s)*(1-t);
    N3=1/8*(1+r)*(1+s)*(1-t);
    N4=1/8*(1-r)*(1+s)*(1-t);
    N5=1/8*(1-r)*(1-s)*(1+t);
    N6=1/8*(1+r)*(1-s)*(1+t);
    N7=1/8*(1+r)*(1+s)*(1+t);
    N8=1/8*(1-r)*(1+s)*(1+t);
end