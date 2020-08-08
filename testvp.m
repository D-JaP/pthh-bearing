clc
clear all
m=100;
x=linspace(0,1,m+1);



luoi=1:(m+1);
luoi_phantu=1:m;
t=zeros(2,m);
for i=1:m
t(1,i)=luoi(i);
t(2,i)=luoi(i+1);
end 
mtA=zeros(m+1);
vtb=zeros(m+1,1);
for i=1:m
    for j=1:3
        switch j
            case 1
                xi=-0.7745966692;
                W(j)=0.555555555;
            case 2
                xi=0;
                W(j)=0.888888888;
            case 3
                xi=0.7745966692;
                W(j)=0.555555555;
        end
        
    L=abs(x(i)-x(i+1));
    A=1/2*[1/L -1/L;-1/L 1/L]*W(j);
    
    mtA(i:i+1,i:i+1)=mtA(i:i+1,i:i+1)+A;
    
    [N1,N2]=shapefunc2D(xi);
    b=-1/3* [N1;N2]*((N1*x(i)+N2*x(i+1))+[N1 N2]*[8;8])*L/2*W(j);
    vtb(i:i+1,1)=vtb(i:i+1,1)+b;
    end
end
vtp=zeros(m+1,1);
vtp(2:m)=mtA(2:m,2:m)\vtb(2:m);

plot(x,vtp)

