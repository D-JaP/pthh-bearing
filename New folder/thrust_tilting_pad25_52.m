clc
clear all
close all
thu=1;
kchia_m=40;kchia_n=40;%so khoang chia
snut=(kchia_m+1)*(kchia_n+1);
spt=kchia_m*kchia_n;

% data=2;
% beta=20/360*2*pi;
% r1=1.110;
% r2=1.685;
% nuy=0.004043;%do n0hot dong luc hoc (N/m^2.s)
% omega=107/60*(2*pi);%toc do goc (rad/s)
% theta_p=11.415/180*pi;
% r_p=1.4175;
% X_p=r_p*sin(theta_p);
% Y_p=r_p*cos(theta_p);
% Wz=1061667;

% data=3;
% beta=20/360*2*pi;
% r1=0.875;
% r2=1.330;
% nuy=0.00565;%do n0hot dong luc hoc (N/m^2.s)
% omega=157/60*(2*pi);%toc do goc (rad/s)
% theta_p=11.652/180*pi;
% r_p=1.1025;
% X_p=r_p*sin(theta_p);
% Y_p=r_p*cos(theta_p);
% Wz=321667;

% Wz=59592;
% beta=55/360*2*pi;
% r1=0.02;
% r2=0.05;
% nuy=0.039;%do nhot dong luc hoc (N/m^2.s)
% omega=3000/60*(2*pi);%toc do goc (rad/s)

% beta=50/360*2*pi;
% r1=0.05715;
% r2=0.1143;
% nuy=0.039;%do nhot dong luc hoc (N/m^2.s)
% omega=1500/60*(2*pi);%toc do goc (rad/s)

data=1;
beta=28/360*2*pi;
r1=0.1875;
r2=0.3225;
nuy=0.0252;%do n0hot dong luc hoc (N/m^2.s)
omega=3000/60*(2*pi);%toc do goc (rad/s)
theta_p=17.38/180*pi;
r_p=0.255;
X_p=r_p*sin(theta_p);
Y_p=r_p*cos(theta_p);
Wz=59592;
% Wz=48766
%% danh so nut
p_so=1;
luoi=zeros(kchia_m+1,kchia_n+1);
for j=1:kchia_n+1
    if mod(j,2)~=0
        for i=1:kchia_m+1
            luoi(i,j)=p_so;
            p_so=p_so+1;
        end
    end
    if mod(j,2)==0
        for i=1:kchia_m+1
            luoi(i,j)=p_so;
            p_so=p_so+1;
        end
    end
end
%% danh so phan tu
p_so2=1;
luoi_phantu=zeros(kchia_m,kchia_n);
for j=1:kchia_n
    %     if mod(i,2)~=0
    for i=1:kchia_m
        luoi_phantu(i,j)=p_so2;
        p_so2=p_so2+1;
    end
    %     end
    %     if mod(i,2)==0
    %         for j=2*n:-1:1
    %             luoi_phantu(j,i)=p_so2;
    %             p_so2=p_so2+1;
    %         end
    %     end
end

%% toa do va luoi phan tu
p=zeros(2,snut);
pp=zeros(2,snut);
for i=1:snut
    [row,col]=find(luoi==i);
    p(1,i)=(col-1)*beta/kchia_n;
    p(2,i)=(r2-(r2-r1)/kchia_m*(row-1));
end

pp(1,:)=p(2,:).*cos(pi/2-p(1,:));
pp(2,:)=p(2,:).*sin(pi/2-p(1,:));
t=zeros(4,spt);
for i=1:kchia_m*kchia_n
    [row,col]=find(luoi_phantu==i);
%     if (col~=n)
        t(1:4,i)=[luoi((row+1),col),luoi((row+1),col+1),luoi(row,col+1),luoi(row,col)];
%     else
%         t(1:4,i)=[luoi((row+1),col),luoi((row+1),1),luoi(row,1),luoi(row,col)];
%     end
end


%%
   
    while(1)
        if (thu==1)
            if (data==1)
                hp=1e-4;
                alpha_r=1e-3;
                alpha_theta=0;
            end
            if (data==2)
                hp=5e-5;
                alpha_r=5e-5;
                alpha_theta=0;
            end
            if (data==3)
                hp=5e-5;
                alpha_r=1e-4;
                alpha_theta=0;
            end
        else

        hp=hp_new;
        alpha_r=alpha_r_new;
        alpha_theta=alpha_theta_new;
            
        end
        
        
        
    dem_ib=1;
    dem_ia=1;
    h=zeros(kchia_m+1,kchia_n+1);
    for i=1:kchia_m+1
        for j=1:kchia_n+1
%             h(i,j)=1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/beta);
%             h(i,j)= (1 + (h1 - h2)*(1 - p(1,luoi(i,j))/beta)/h2);
%             h(i,j)=(1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/sqrt(pp(1,luoi(i,j))^2 + pp(2,luoi(i,j))^2)));
            h(i,j)=(hp+p(2,luoi(i,j))*sin(theta_p-p(1,luoi(i,j)))*sin(alpha_r)+(r_p-p(2,luoi(i,j))*cos(theta_p-p(1,luoi(i,j))))*sin(alpha_theta))/hp;
            
            if (i==1||j==1||i==kchia_m+1||j==kchia_n+1)
            I_boundry(dem_ib)=luoi(i,j);
            dem_ib=dem_ib+1;
            else
            I_active(dem_ia)=luoi(i,j);
            dem_ia=dem_ia+1;
            end
        end
    end
    if (sum(h(:,:)<0)~=0)
        disp('error h')
    end
%%    
    H(luoi(:,:))=h(:,:);% chuyen array 2d h thanh H 1d
     if (sum(H(:)<0)~=0)
        disp('error h')
    end
    % tich phan so
    mtA=zeros(snut,snut);
    A_fix=zeros(snut,snut);
    mt_dA_dhp=zeros(snut,snut);% Partial derivatives A/hp
    mt_dA_dalpha_r=zeros(snut,snut);% Partial derivatives A/alpha_r
    mt_dA_dalpha_theta=zeros(snut,snut);% Partial derivatives A/alpha_theta
    
    b=zeros(snut,1);
    bh=zeros(snut,1);
    dB_dhp=zeros(snut,1);% Partial derivatives B/hp
    dB_dalpha_r=zeros(snut,1);% Partial derivatives B/alpha_r
    dB_dalpha_theta=zeros(snut,1);% Partial derivatives B/alpha_theta
    
    unit_vt=zeros(snut,1);%unit vecto for Gauss integral
    
    vt_X=zeros(snut,1);%unit vecto X for Gauss integral
    vt_Y=zeros(snut,1);%unit vecto Y for Gauss integral
    %tich phan so Gauss
    
    for i=1:spt
            theta=p(1,t(1:4,i));
            r=p(2,t(1:4,i));
            xP=pp(1,t(1:4,i));
            yP=pp(2,t(1:4,i));
        for j=1:9
            [xi,eta,W]=tpso_Gauss_2d(3,j);
            [N1,N2,N3,N4]=hamdang(xi,eta);
            N=[N1 ;N2 ;N3 ;N4];
            dNij=[ -1/4*(1-eta) -1/4*(1-xi);...
                1/4*(1-eta) -1/4*(1+xi);...
                1/4*(1+eta)  1/4*(1+xi);...
                -1/4*(1+eta)  1/4*(1-xi)];
            
            J=Jctugiac(theta,r,eta,xi);
            Jn=inv(J);
            
            dNx(1:4)=Jn(1,1)*dNij(1:4,1)+Jn(1,2)*dNij(1:4,2);
            dNz(1:4)=Jn(2,1)*dNij(1:4,1)+Jn(2,2)*dNij(1:4,2);
            
            A=N'*theta';%bien
            B=N'*r';%bien
            C=N'*xP';%toa do de cac
            D=N'*yP';%toa do de cac
            
%             HH3=(N'*H(t(1:4,i))')^3;
%             HH2=(N'*H(t(1:4,i))')^2;
%             HH=(N'*H(t(1:4,i))');
            THETA_p=(N'*[theta_p;theta_p;theta_p;theta_p]);%%%

            HH=N'*(H(t(1:4,i)))';
            HH2=HH^2;
            HH3=HH^3;
            %tinh cac ma tran do cung
            A_i=HH3*(1/B.*dNx'*dNx+B*(dNz'*dNz))*B*abs(det(J))*W(j);%toa do cuc
            dA_dhp_i=3*HH2*(B*sin(-theta_p + A)*sin(alpha_r) + sin(alpha_theta)*cos(-theta_p + A)*B - sin(alpha_theta)*r_p)/hp^2*(1/B.*dNx'*dNx+B*(dNz'*dNz))*B*abs(det(J))*W(j);
            dA_dalpha_r_i=3*HH2*-B*sin(-theta_p + A)*cos(alpha_r)/hp*(1/B.*dNx'*dNx+B*(dNz'*dNz))*B*abs(det(J))*W(j);
            dA_dalpha_theta_i=3*HH2*(r_p - B*cos(-theta_p + A))*cos(alpha_theta)/hp*(1/B.*dNx'*dNx+B*(dNz'*dNz))*B*abs(det(J))*W(j);
%             Ah=3*HH2*-(beta - acot(B/A))*h1/(h2^2*beta)*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
%             Ah=3*HH2*-(beta - acot(B/A))*h1/(h2^2*beta)*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j); 
            mtA(t(1:4,i),t(1:4,i))=mtA(t(1:4,i),t(1:4,i))+A_i;
            mt_dA_dhp(t(1:4,i),t(1:4,i))=mt_dA_dhp(t(1:4,i),t(1:4,i))+dA_dhp_i;
            mt_dA_dalpha_r(t(1:4,i),t(1:4,i))=mt_dA_dalpha_r(t(1:4,i),t(1:4,i))+dA_dalpha_r_i;
            mt_dA_dalpha_theta(t(1:4,i),t(1:4,i))=mt_dA_dalpha_theta(t(1:4,i),t(1:4,i))+dA_dalpha_theta_i;
%             mtAh(t(1:4,i),t(1:4,i))=mtAh(t(1:4,i),t(1:4,i))+Ah;
            %tinh vecto tai
%             A=xP;
%             B=yP;
%             dH=h2*delta*(B./(A.^2.*(B.^2./A.^2 + 1).*sqrt(A.^2 + B.^2)) + acot(B./A).*A./(A.^2 + B.^2).^(3/2));
            
            
            
            
%             dH=(h2-h1)/h2*B/(A^2*(B^2/A^2 + 1)*beta);
            
            %dH=delta*(-B/(A^2*(B^2/A^2 + 1)*sqrt(A^2 + B^2)) + acot(B/A)*A/(A^2 + B^2)^(3/2));
            
            
            %Bj=6*N*(N'*dH')*(N'*sqrt(A.^2+B.^2)')*abs(det(J))*W(j);
%            Bj=-6*nuy*omega/h2^2*N*dH*r*abs(det(J))*W(j);
%            Bj=-6*nuy*omega/h2^2*N*-(h1 - h2)*B^2/((A^2 + B^2)*h2*beta)*abs(det(J))*W(j);
            %Bj=-6*nuy*omega/h2^2*N*r*-(h1 - h2)*B/(h2*(A^2 + B^2)*beta)*abs(det(J))*W(j);
%             B_i=-6*nuy*omega*N*B*( -B*cos(theta_p-A)*sin(alpha_r) +(r_p-B*sin(theta_p-A))*sin(alpha_theta) )*B*abs(det(J))*W(j)/hp^3;
            B_i=-6*nuy*omega*B*N*(dNx*H(t(1:4,i))')/hp^2*abs(det(J))*B*W(j);
            dB_dhp_i=-6*nuy*omega*N*B*3*B*(sin(alpha_r)*cos(-theta_p + A) - sin(alpha_theta)*sin(-theta_p + A))/hp^4*B*abs(det(J))*W(j);
            dB_dalpha_r_i=-6*nuy*omega*N*B*-B*cos(-theta_p + A)*cos(alpha_r)/hp^3*B*abs(det(J))*W(j);
            dB_dalpha_theta_i=-6*nuy*omega*N*B*B*sin(-theta_p + A)*cos(alpha_theta)/hp^3*B*abs(det(J))*W(j);
            %bhi=-6*nuy*omega*N*(A/sqrt(A^2 + B^2))*abs(det(J))*W(j)/h2^2;
%             bhi=-6*nuy*omega*N*B*(3*h1 - 2*h2)/(sqrt(A^2 + B^2)*h2^4*beta)*abs(det(J))*W(j);
            
            b(t(1:4,i))=b(t(1:4,i))+B_i;
            dB_dhp(t(1:4,i))=dB_dhp(t(1:4,i))+dB_dhp_i;
            dB_dalpha_r(t(1:4,i))=dB_dalpha_r(t(1:4,i))+dB_dalpha_r_i;
            dB_dalpha_theta(t(1:4,i))=dB_dalpha_theta(t(1:4,i))+dB_dalpha_theta_i;
            
%             X=N*C*B*abs(det(J))*W(j);
%             Y=N*D*B*abs(det(J))*W(j);
            X=N*B^2*sin(A-theta_p)*abs(det(J))*W(j);
            Y=N*B*(B*cos(A-theta_p)-r_p)*abs(det(J))*W(j);
            
            vt_X(t(1:4,i))=vt_X(t(1:4,i)) + X;%vecto x cho chuyen doi nodal forces
            vt_Y(t(1:4,i))=vt_Y(t(1:4,i)) + Y;%vecto y cho chuyen doi nodal forces
            
            unit_vti=N*B*abs(det(J))*W(j);
            unit_vt(t(1:4,i))=unit_vt(t(1:4,i)) + unit_vti; % tao vecto don vi cho tich phan gauss
            
        end
    end
for i=1:snut
    if (sum(I_active==i)==0) % i khong thuoc I active
        A_fix(i,i)=1;
    end
end

mt_dA_dhp_fix=zeros(snut,snut);
mt_dA_dalpha_r_fix=zeros(snut,snut);
mt_dA_dalpha_theta_fix=zeros(snut,snut);
% dp_dhp=zeros(snut,1);
% dp_dalpha_r=zeros(snut,1);
% dp_dalpha_theta=zeros(snut,1);

A_fix(I_active,I_active)=mtA(I_active,I_active);clear mtA
mt_dA_dhp_fix(I_active,I_active)=mt_dA_dhp(I_active,I_active);clear mt_dA_dhp
mt_dA_dalpha_r_fix(I_active,I_active)=mt_dA_dalpha_r(I_active,I_active);clear mt_dA_dalpha_r
mt_dA_dalpha_theta_fix(I_active,I_active)=mt_dA_dalpha_theta(I_active,I_active);clear mt_dA_dalpha_theta

mt_dA_dhp_fix(I_active,I_boundry)=0;mt_dA_dhp_fix(I_boundry,I_active)=0;
mt_dA_dalpha_r_fix(I_active,I_boundry)=0;mt_dA_dalpha_r_fix(I_boundry,I_active)=0;
mt_dA_dalpha_theta_fix(I_active,I_boundry)=0;mt_dA_dalpha_theta_fix(I_boundry,I_active)=0;

b(I_boundry)=0;%still not fix but straightly to b
dB_dhp(I_boundry)=0;
dB_dalpha_r(I_boundry)=0;
dB_dalpha_theta(I_boundry)=0;

% bh(I_boundry)=0;
vtp=(A_fix\b)';% tinh vecto ap suat

Fz=unit_vt'*vtp';

M_X=vt_X'*vtp';
M_Y=vt_Y'*vtp';

% XcF=M_X/Fz
% YcF=M_Y/Fz

dp_dhp=(A_fix\(dB_dhp-mt_dA_dhp_fix*vtp'));
dp_dalpha_r=(A_fix\(dB_dalpha_r-mt_dA_dalpha_r_fix*vtp'));
dp_dalpha_theta=(A_fix\(dB_dalpha_theta-mt_dA_dalpha_theta_fix*vtp'));

clear A_fix dB_dhp dB_dalpha_r dB_dalpha_theta mt_dA_dhp_fix mt_dA_dalpha_r_fix mt_dA_dalpha_theta_fix
dFz_dhp=unit_vt'*dp_dhp; 
dX_dhp=vt_X'*dp_dhp;
dY_dhp=vt_Y'*dp_dhp;

dFz_dapha_r=unit_vt'*dp_dalpha_r;
dX_dapha_r=vt_X'*dp_dalpha_r;
dY_dapha_r=vt_Y'*dp_dalpha_r;

dFz_dapha_theta=unit_vt'*dp_dalpha_theta;
dX_dapha_theta=vt_X'*dp_dalpha_theta;
dY_dapha_theta=vt_Y'*dp_dalpha_theta;

clear unit_vt vt_X vt_Y dp_dhp dp_dalpha_r dp_dalpha_theta

D=[dFz_dhp dFz_dapha_r dFz_dapha_theta;dX_dhp dX_dapha_r dX_dapha_theta;dY_dhp dY_dapha_r dY_dapha_theta];
% D=[ dX_dapha_r dX_dapha_theta; dY_dapha_r dY_dapha_theta];
new_variables=[hp;alpha_r;alpha_theta]-D\[Fz-Wz;M_X-0;M_Y-0]
% D=[dFz_dhp dFz_dapha_theta ;dX_dhp dX_dapha_theta ;dY_dhp dY_dapha_theta ];
% new_variables=[hp;alpha_theta]-D\[Fz-Wz;XcF-X_p;YcF-Y_p]
hp_new=new_variables(1);
alpha_r_new=new_variables(2);
r_table(thu)=alpha_r;
alpha_theta_new=new_variables(3);
theta_table(thu)=alpha_theta;
% (Fz-Wz)/Fz
disp('max pressure')

P_table(thu)=max(vtp);

disp('min film thickness')
min(H)*hp

disp('Fz Mx My')
Fz
M_X
M_Y

F_table(thu)=Fz;
X_table(thu)=M_X;
Y_table(thu)=M_Y;
Hp_table(thu)=hp;

if (    norm(alpha_r_new-alpha_r)/norm(alpha_r_new)<=1e-5&&...
        norm(alpha_theta_new-alpha_theta)/norm(alpha_theta_new)<=1e-5&&...
        norm(Fz-Wz)/norm(Fz)<=1e-5 &&...
        abs(M_X)<=1e-5&&...
        abs(M_Y)<=1e-5)
    
    disp('done')
break
end
thu=thu+1;
    end
%%
    % cd=linspace(0,1,m+1);
    % mtp(:,:)=vtp(luoi(:,:));
    % mtq(:,:)=vtq(luoi(:,:));
    % [X,Y]=meshgrid(cd,phi);
    % figure(1)
    % surf(X',Y',mtp)
    % figure(2)
    % surf(X',Y',mtq)
    
    
    
    % vtp_giua=vtp(luoi(m/2,:));
    % figure(3)
    % plot(phi,vtp_giua)
    % phi(find(vtp_giua==max(vtp_giua)));

    % for i=1:snut
    %     if (sum(I_active==i)==0) % i khong thuoc I active
    %         Am_X(i,i)=1;
    %         Am_Y(i,i)=1;
    %     end
    % end
    
    %tinh dp/dx, dp/dy
%     pxy=Am\[-Am_X*vtp+Bm_X -Am_Y*vtp+Bm_Y];
%     
%     pnx=pxy(:,1);
%     pny=pxy(:,2);
%     
%     D_jc=-[vt_S'*pnx vt_S'*pny;...
%         vt_R'*pnx vt_R'*pny];
%     
%     fx=-vt_S'*vtp
%     fy=-vt_R'*vtp
%     
%     u_new=u-D_jc\([fx;fy]-[Wx;Wy])
%     
%     if ((norm(u_new-u)/norm(u_new))<=1e-5&&(norm([fx;fy]-[Wx;Wy])/norm([Wx;Wy]))<1e-5)
%         disp('done')
%         break
%     end
%     
%     u=u_new;
%     thu=thu+1;
%%


mtp=zeros(kchia_m+1,kchia_n+1);

mtp(:,:)=vtp(luoi(:,:));

% [X,Y,Z]=meshgrid(pp(1,:),pp(2,:),vtp)
% 
% s=surf(X,Y,Z);
% figure(2)
% plot(pp(1,:),pp(2,:),'+')
figure(1)
x=pp(1,:)';y=pp(2,:)';z=vtp';
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ; 
yi = dt.Points(:,2) ; 
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
trisurf(tri,xi,yi,zi) 
daspect([1 1 max(abs(vtp))])
figure(2)
x=pp(1,:)';y=pp(2,:)';z=H';
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ; 
yi = dt.Points(:,2) ; 
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
trisurf(tri,xi,yi,zi) 
daspect([1 1 max(abs(H))])
%%
figure(3)
[row,col]=find(mtp==max(vtp));
kc=(0:kchia_n)*beta/kchia_n/pi*180;
kc1=mtp(kchia_m/2+1,1:kchia_n+1)*1e-6;
kc2=mtp(kchia_m*0.8+1,1:kchia_n+1)*1e-6;
kc3=mtp(kchia_m*0.2+1,1:kchia_n+1)*1e-6;

plot(kc,kc1,'-',kc,kc2,'--',kc,kc3,'-.')
c1=imregionalmax(mtp(kchia_m/2+1,1:kchia_n+1)*1e-6);
c2=imregionalmax(mtp(kchia_m*0.8+1,1:kchia_n+1)*1e-6);
c3=imregionalmax(mtp(kchia_m*0.2+1,1:kchia_n+1)*1e-6);
hold on
plot(kc(c1),kc1(c1),'*',kc(c2),kc2(c2),'*',kc(c3),kc3(c3),'*')
kc1(c1)
kc2(c2)
kc3(c3)
xlabel('Circumferential Direction (deg)','fontweight','bold','fontsize',11,'color','black');
ylabel('Pressure - P(MPa)','fontweight','bold','fontsize',11,'color','black')
grid on

if (data==1)
    axis([0 28 0 9])
    xticks([0:2:28])
    legend('r=255 m','r=295.5 m','r=214.5 m')
end 
if (data==2)
    axis([0 28 0 9])
    xticks([0:2:28])
    legend('r=1.3975 m','r=2.1580 m','r=1.8070 m')
end

if (data==3)
    axis([0 20 0 4.5])
    xticks([0:2:20])
    legend('r=1.1025 m','r=1.2390 m','r=0.9660 m')
end
%%
% figure(4)
% 
% plot(kchia_n:-1:0,mtp(1:kchia_m+1,kchia_n*0.7+1),kchia_n:-1:0,mtp(1:kchia_m+1,kchia_n*0.6+1),kchia_n:-1:0,mtp(1:kchia_m+1,kchia_n*0.1+1))
% xlabel('Circumferential Direction(deg)','fontweight','bold','fontsize',11,'color','black');
% ylabel('Pressure - P(Pa)','fontweight','bold','fontsize',11,'color','black')
% grid on
% axis([0 28 0 9e6])
% xticks([0:2:28])
%%
if (data==1)
figure(10)
plot((0:kchia_n)*beta/kchia_n/pi*180,1e-6*mtp(kchia_m/2+1,1:kchia_n+1),'-k')
x=[0    0.921316    1.88269    2.80401    3.7253    4.6867    5.5679    6.5293    7.4506    8.4120    9.3333   10.2546   11.2160   12.1774   13.0937   14.0200   14.9413   15.8627   16.8240   17.7454   18.6667   19.6230   20.5494   21.4707   22.3920   23.3534   24.2747   25.2361   26.1574   27.0787 28];
y=[0    0.404816    0.82653    1.248980    1.6714   2.0934    2.4796    2.9020    3.3245    3.7286    4.1510    4.5913    5.0327    5.4735    5.9142    6.337    6.7776    7.1816    7.5857    7.9347   8.2469    8.4857    8.6510    8.6878    8.5592    8.2286    7.6224    6.6306    5.1612    3.0306 0];
hold on
plot(x,y,'-or')
xlabel('Circumferential Direction(deg)','fontweight','bold','fontsize',11,'color','black');
ylabel('Pressure - P(Pa)','fontweight','bold','fontsize',11,'color','black')
grid on
axis([0 28 0 9])
legend('Present Study','[4]')
xticks([0:2:28])
end 
%%
figure(6)
x=pp(1,:)';y=pp(2,:)';z=vtp';
for i=1:kchia_m
line(x(luoi(i,:)),y(luoi(i,:)),z(luoi(i,:)))
end
for i=1:kchia_n
line(x(luoi(:,i)),y(luoi(:,i)),z(luoi(:,i)))
end
a=[luoi(1,:)';luoi(:,kchia_n+1)];
line(x(a),y(a),z(a))
b=[luoi(:,1);luoi(kchia_m+1,:)'];
line(x(b),y(b),z(b))

axis square
%%
figure(7)
theta=linspace(pi/2,pi/2-beta,kchia_n+1);
r=linspace(r1,r2,kchia_m+1);
[theta,r]=meshgrid(theta,r);
[x,y,z]=pol2cart(theta,r,mtp*1e-6);
surf(x,y,z)
c=imregionalmax(z);
hold on 
k1=plot3(x(c),y(c),z(c),'xy');
k2=plot3(X_p,Y_p,z(c),'+k');
legend([k1,k2],'Max pressure','Pivot')

hold off

set(gca,'DataAspectRatio',[1 1 5])

% axis tight
xlabel('x (m)','fontweight','bold','fontsize',11,'color','black');
ylabel('y (m)','fontweight','bold','fontsize',11,'color','black');
zlabel('Pressure - P(Pa)','fontweight','bold','fontsize',11,'color','black')
alpha 0.8
colormap(jet)
% i = rotate3d;set(i,'ActionPreCallback',...
%     'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(i,'ActionPostCallback',...
%     'set(gcf,''windowbuttonmotionfcn'','''')')
% lighting phong
%% 
figure(8)
contour3(x,y,z/1e6,[  2  4  6 8],'ShowText','on')
set(gca,'DataAspectRatio',[1 1 2])
hold on

a=[luoi(1,:)';luoi(:,kchia_n+1)];
line(x(a),y(a),z(a)*0)
b=[luoi(:,1);luoi(kchia_m+1,:)'];
line(x(b),y(b),z(b)*0)

axis equal
axis([-0.1 0.6 1 1.75  0 10])
% plot3()
%%
figure(9)
theta=linspace(pi/2,pi/2-beta,kchia_n+1);
r=linspace(r1,r2,kchia_m+1);
[theta,r]=meshgrid(theta,r);
[x,y,z]=pol2cart(theta,r,h*hp*1e6);
surf(x,y,z)
set(gca,'DataAspectRatio',[1 1 1e7*2])
% axis tight
hold on 

k2=plot3(X_p,Y_p,z(c),'+k');
legend(k2,'Pivot')

hold off
xlabel('x (m)','fontweight','bold','fontsize',11,'color','black');
ylabel('y (m)','fontweight','bold','fontsize',11,'color','black');
zlabel('Pessure - P(Pa)','fontweight','bold','fontsize',11,'color','black')
alpha 0.8
colormap(flipud(jet))
%% 
% figure(1)
% s=surf(X',Y',mtp*6*nuy*omega*(R/c)^2);
% axis tight
% xlabel('Bearing length - L(mm)','fontweight','bold','fontsize',11,'color','black');
% ylabel('Cicumference - \theta (rad)','fontweight','bold','fontsize',11,'color','black');
% zlabel('Pressure - P(Pa)','fontweight','bold','fontsize',11,'color','black')
% yticks([0:60:360])
% xticks([0:10:50])
% view(45,45)
% alpha 0.8
% pbaspect([1 1.5 1])
% i = rotate3d;set(i,'ActionPreCallback',...
%     'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(i,'ActionPostCallback',...
%     'set(gcf,''windowbuttonmotionfcn'','''')')
% lighting phong
% 
% set(s,'edgecolor',[0 0 0.4],'meshstyle','both','linewidth',.15);
% figure(2)
% surf(X',Y',mtq)
% figure(3)
% plot(phi/pi*180,h)
% 
% figure(4)
% mtR(:,1:n)=vt_R(luoi(:,1:n));
% mtS(:,1:n)=vt_S(luoi(:,1:n));
% mtR(:,n+1)=mtR(:,1);
% mtS(:,n+1)=mtS(:,1);
% surf(X',Y',mtR)
% figure(5)
% surf(X',Y',mtS)
%%
% figure(6)
% o=1+(x*cos(phi)+y*sin(phi));
% yyaxis left
% plot(phi(1:n+1)/pi*180,mtp(m/2+1,1:n+1)*6*nuy*omega*(R/c)^2,'k-');
% hold on
% xticks([0:60:360])
%
% xlim([0,360])
% plot(phi(1:n+1)/pi*180,mtp(11,1:n+1)*6*nuy*omega*(R/c)^2,'k-');
% plot(phi(1:n+1)/pi*180,mtp(5,1:n+1)*6*nuy*omega*(R/c)^2,'k-');
% hold off
% xlabel('cicumference - \theta (^{o})','fontweight','bold','fontsize',11,'color','black');
% ylabel('pressure - P(Pa)','fontw $x^2+e^{\pi i}$ eight','bold','fontsize',11,'color','black')
% ylim([0 100000])
% grid on 
% yyaxis right
% plot(phi(1:n+1)/pi*180,o*c*1000,'k-');
% ylabel('film thickness  - h(mm)','fontweight','bold','fontsize',11,'color','black')
% ylim([0,0.07])
%%
figure(11)
theta=linspace(0,beta,kchia_n+1);
[row,col]=find(mtp==max(vtp));
plot(theta,h(row,:)*hp*1e6,'-r',theta,h(kchia_m/2+1,:)*hp*1e6,'--')
hold on
legend('P max','P mid')

figure(12)
theta=linspace(0,beta,kchia_n+1);
[row,col]=find(mtp==max(vtp));
plot(theta,mtp(row,:),'-r',theta,mtp(kchia_m/2+1,:),'--')
hold on
legend('r=','r=')
%%
figure(13)
theta=linspace(pi/2,pi/2-beta,kchia_n+1);
r=linspace(r1,r2,kchia_m+1);
[theta,r]=meshgrid(theta,r);
[x,y,z]=pol2cart(theta,r,mtp*1e-6);
surf(x,y,z)


set(gca,'DataAspectRatio',[1 1 5])

% axis tight
xlabel('x (m)','fontweight','bold','fontsize',11,'color','black');
ylabel('y (m)','fontweight','bold','fontsize',11,'color','black');
zlabel('Pressure - P(MPa)','fontweight','bold','fontsize',11,'color','black')
alpha 0.8
colormap(jet)
% i = rotate3d;set(i,'ActionPreCallback',...
%     'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(i,'ActionPostCallback',...
%     'set(gcf,''windowbuttonmotionfcn'','''')')
% lighting phong
%%
table=[Hp_table' r_table' theta_table' P_table' F_table' X_table' Y_table']

[row,col]=find(mtp==max(vtp));
disp('r at max')
(kchia_m-(row-1))/kchia_m*(r2-r1)+r1
disp('theta at max')
((col-1))/kchia_n*beta/pi*180


