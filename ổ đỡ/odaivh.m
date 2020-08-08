clc
clear all
m=100;
phi=linspace(0,2*pi,m+1);
epsilon=0.3;
h=1+epsilon*cos(phi);

luoi=1:(m);
luoi_phantu=1:m;
t=zeros(2,m);
for i=1:m
    if i~=m
        t(1,i)=luoi(i);
        t(2,i)=luoi(i+1);
    end
    
    if i==m
        t(1,i)=luoi(i);
        t(2,i)=luoi(1);
    end
end 
mtA=zeros(m);
vtb=zeros(m,1);
vt_R=zeros(m,1);
vt_S=zeros(m,1);
L=abs(phi(1)-phi(1+1));
for i=1:m
    for j=1:2
        [xi,W]=tpso_Gauss_1d(2,j);
        [N1,N2]=shapefunc2D(xi);
        N=[N1;N2];
%     L=abs(phi(i)-phi(i+1));
    dNx=inv(L/2)*
    A=1/2*[1/L -1/L;-1/L 1/L]*W(j)*([N1 N2]*h(t(:,i))').^3;
    
    mtA(t(:,i),t(:,i))=mtA(t(:,i),t(:,i))+A;
    if i~=m
    b=6*[N1;N2]*epsilon*sin(N1*phi(t(1,i))+N2*phi(t(2,i)))*L/2*W(j);
    else
    b=6*[N1;N2]*epsilon*sin(N1*phi(t(1,i))+N2*phi(m+1))*L/2*W(j);
    end
    vtb(t(:,i),1)=vtb(t(:,i),1)+b;
    if i~=m
    RR=N*cos(N'*phi(t(:,i))')*L/2*W(j);
    SS=N*sin(N'*phi(t(:,i))')*L/2*W(j);
    else
    RR=N*cos(N'*[phi(t(1,i));phi(m+1)])*L/2*W(j);
    SS=N*sin(N'*[phi(t(1,i));phi(m+1)])*L/2*W(j);
    end
    vt_R(t(:,i),1)=vt_R(t(:,i),1)+RR;
    vt_S(t(:,i),1)=vt_S(t(:,i),1)+SS;
    end
end

vtp=zeros(m,1);
[fi_cav,gamma_cav]=tinh_fi_cav(epsilon);
for i=length(phi):-1:1
 
   if (phi(i)<=fi_cav)
       vtp(i:length(vtp))=0;
       mesh_point=i;
       break
   end
   
end

vtp(2:mesh_point)=mtA(2:mesh_point,2:mesh_point)\vtb(2:mesh_point);


disp('P max :') ,max(vtp)
disp('fi max:') ,phi(find(vtp==max(vtp)))
disp('Wx') ,Wx=(-3*epsilon)/(1-epsilon^2)*(1-cos(gamma_cav))^2/(1-epsilon*cos(gamma_cav))
Wy=-6/sqrt(1-epsilon^2)*(gamma_cav*cos(gamma_cav)-sin(gamma_cav))/(1-epsilon*cos(gamma_cav))
disp('Wx Wy:')
-vt_R'*vtp
-vt_S'*vtp
plot(phi/pi*180,h)
% figure(2)
% plot(phi(1:m+1)/pi*180,vt_R*0.035)
figure(3)
dp_dphi=(vtp(2:mesh_point)-vtp((2:mesh_point)-1))./(2*pi/m)



L=abs(phi(1)-phi(1+1));
dp_dphi2=zero(m,1);
for i=1:m
    for j=1:2
        [xi,W]=tpso_Gauss_1d(2,j);
        [N1,N2]=shapefunc2D(xi);
        N=[N1;N2];
%     L=abs(phi(i)-phi(i+1));
    dp_dphi2_i=N*
    dp_dphi2(t(:,i),1)=dp_dphi2(t(:,i),1)+dp_dphi2_i;
    
    end
end