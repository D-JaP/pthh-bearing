clc
clear all
close all


thu=1;
input_ochan_tilting_nhiet
%%
global snut snut3D spt2D spt3D 
snut=(kchia_m+1)*(kchia_n+1);
snut3D=(kchia_m+1)*(kchia_n+1)*(kchia_q+1);
spt2D=kchia_m*kchia_n;
spt3D=kchia_m*kchia_n*kchia_q;
%% danh so nut luoi 2D
global luoi
[luoi]=danhsonutluoi2D(kchia_m,kchia_n);
%% danh so nut luoi 3D va 3Dto2D
global luoi3D
[luoi3D]=danhsonutluoi3D(kchia_m,kchia_n,kchia_q);% danh so nut


%% danh so phan tu 2D 3D
global luoi_phantu luoi_phantu3D luoi_phantu3D_all
luoi_phantu=luoiphantu2D(kchia_m,kchia_n);

luoi_phantu3D=luoiphantu3D(kchia_m,kchia_n,kchia_q);
luoi_phantu3D_all=luoiphantu3D(kchia_m,kchia_n,kchia_q*2);
%% toa do va luoi phan tu 2D

p=zeros(2,snut);
pp=zeros(2,snut);
for i=1:snut
    [row,col]=find(luoi==i);
    p(1,i)=(col-1)*beta/kchia_n;
    p(2,i)=(r2-(r2-r1)/kchia_m*(row-1));
end

pp(1,:)=p(2,:).*cos(pi/2-p(1,:));
pp(2,:)=p(2,:).*sin(pi/2-p(1,:));


%% chuyen cac nut vao phan tu
global t
t=chuyennutvaophantu(spt2D,luoi_phantu,kchia_m,kchia_n,luoi);


%%
luoi_pad3D=luoi3D+(kchia_m+1)*(kchia_n+1)*kchia_q;
snut_all=max(luoi_pad3D,[],'all');
luoi_all=zeros(kchia_m+1,kchia_n+1,kchia_q*2+1);
luoi_all(:,:,1:kchia_q+1)=luoi3D(:,:,:);
luoi_all(:,:,(kchia_q+1):(kchia_q*2+1))=luoi_pad3D(:,:,:);
    %%
%     global h1 h2 deltah h2_new delta
%     while(1)
        if (thu==1)

        hp=100e-6;
        alpha_r=0.0003;
        alpha_theta=0;
        else

        hp=100e-6;
        alpha_r=alpha_r_new;
        alpha_theta=alpha_theta_new;
            
        end
        
        
%% tinh h
dem_ib=1;
dem_ia=1;
h=zeros(kchia_m+1,kchia_n+1);
%     phuongchieucao=zeros(kchia_n+1,kchia_n+1,kchia_n+1);%%%%%%%%%%%%%%%%%%
%     PHUONGCHIEUCAO=zeros((kchia_n+1)*(kchia_n+1),kchia_n+1);%%%%%%%%%%%%%%%%%%
    for i=1:kchia_m+1
        for j=1:kchia_n+1
%             h(i,j)=1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/beta);
%             h(i,j)= (1 + (h1 - h2)*(1 - p(1,luoi(i,j))/beta)/h2);
%             h(i,j)=(1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/sqrt(pp(1,luoi(i,j))^2 + pp(2,luoi(i,j))^2)));
            h(i,j)=(hp+p(2,luoi(i,j))*sin(theta_p-p(1,luoi(i,j)))*sin(alpha_r)+(r_p-p(2,luoi(i,j))*cos(theta_p-p(1,luoi(i,j))))*sin(alpha_theta))/hp;
            dh_dtheta(i,j)=-p(2,luoi(i,j))*cos(-theta_p + p(1,luoi(i,j)))*sin(alpha_r)/hp + p(2,luoi(i,j))*sin(-theta_p + p(1,luoi(i,j)))*sin(alpha_theta)/hp;
            dh_dr(i,j)=-sin(-theta_p + p(1,luoi(i,j)))*sin(alpha_r)/hp - cos(-theta_p + p(1,luoi(i,j)))*sin(alpha_theta)/hp;

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
    Y3D=zeros((kchia_m+1)*(kchia_n+1)*(kchia_q+1),1);
    H3D_global=zeros((kchia_m+1)*(kchia_n+1)*(kchia_q+1),1);
    
    for i=1:kchia_m+1
        for j=1:kchia_n+1
            for k=1:kchia_q+1
                Y3D(luoi3D(i,j,k))=h(i,j)/(kchia_q)*(k-1);
                H3D_global(luoi3D(i,j,k))=h(i,j);
%                 Y3D(luoi3D(i,j,k))/H3D_global(luoi3D(i,j,k))
%                 if k==kchia_q+1
%                     H3D(luoi3D(i,j,k))/H3D_global(luoi3D(i,j,k))
%                 end
            end
        end
    end
    
%%  
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
    dh_dtheta2=zeros(snut,1);
    %tich phan so Gauss

%% tich phan so Gauss
    for i=1:spt2D
            theta=p(1,t(1:4,i));
            r=p(2,t(1:4,i));
            xP=pp(1,t(1:4,i));
            yP=pp(2,t(1:4,i));
        for j=1:9
            [xi,eta,trong_so]=tpso_Gauss_2d(3,j);
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
            
            
            HH=N'*(H(t(1:4,i)))';
            HH2=HH^2;
            HH3=HH^3;
            %tinh cac ma tran do cung
            A_i=                               (HH3*(1/B.*dNx'*dNx+B*(dNz'*dNz)))*B*abs(det(J))*trong_so(j);%toa do cuc
            
            dh_dtheta_i=(dNx*(H(t(1:4,i)))')*B*trong_so(j)*abs(det(J));
            dh_dtheta2(t(1:4,i))=dh_dtheta2(t(1:4,i))+dh_dtheta_i;
            
            dA_dalpha_r_i=               -B*sin(-theta_p + A)*cos(alpha_r)/hp*(3*HH2*(1/B.*dNx'*dNx+B*(dNz'*dNz)))*B*abs(det(J))*trong_so(j);
            dA_dalpha_theta_i=(r_p - B*cos(-theta_p + A))*cos(alpha_theta)/hp*(3*HH2*(1/B.*dNx'*dNx+B*(dNz'*dNz)))*B*abs(det(J))*trong_so(j);
%             Ah=3*HH2*-(beta - acot(B/A))*h1/(h2^2*beta)*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
%             Ah=3*HH2*-(beta - acot(B/A))*h1/(h2^2*beta)*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j); 
            mtA(t(1:4,i),t(1:4,i))=mtA(t(1:4,i),t(1:4,i))+A_i;
            
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
%            Bj=-6*nuy0*omega/h2^2*N*dH*r*abs(det(J))*W(j);
%            Bj=-6*nuy0*omega/h2^2*N*-(h1 - h2)*B^2/((A^2 + B^2)*h2*beta)*abs(det(J))*W(j);
            %Bj=-6*nuy0*omega/h2^2*N*r*-(h1 - h2)*B/(h2*(A^2 + B^2)*beta)*abs(det(J))*W(j);
            
            B_i=-6*nuy0*omega*N*B^2*(-sin(alpha_r)*cos(-theta_p + A) + sin(alpha_theta)*sin(-theta_p + A))/hp*B*abs(det(J))*trong_so(j);
%             B_i=-6*nuy0*omega*B*N*(dNx*H(t(1:4,i))')*abs(det(J))*trong_so(j);
            
            dB_dalpha_r_i=-6*nuy0*omega*N*-B^2*cos(-theta_p + A)*cos(alpha_r)/hp*B*abs(det(J))*trong_so(j);
            dB_dalpha_theta_i=-6*nuy0*omega*N*B^2*sin(-theta_p + A)*cos(alpha_theta)/hp*B*abs(det(J))*trong_so(j);
            %bhi=-6*nuy0*omega*N*(A/sqrt(A^2 + B^2))*abs(det(J))*W(j)/h2^2;
%             bhi=-6*nuy0*omega*N*B*(3*h1 - 2*h2)/(sqrt(A^2 + B^2)*h2^4*beta)*abs(det(J))*W(j);
            
            b(t(1:4,i))=b(t(1:4,i))+B_i;
            
            dB_dalpha_r(t(1:4,i))=dB_dalpha_r(t(1:4,i))+dB_dalpha_r_i;
            dB_dalpha_theta(t(1:4,i))=dB_dalpha_theta(t(1:4,i))+dB_dalpha_theta_i;
            
            X=N*C*B*abs(det(J))*trong_so(j);
            Y=N*D*B*abs(det(J))*trong_so(j);
            
            vt_X(t(1:4,i))=vt_X(t(1:4,i)) + X;%vecto x cho chuyen doi nodal forces
            vt_Y(t(1:4,i))=vt_Y(t(1:4,i)) + Y;%vecto y cho chuyen doi nodal forces
            
            unit_vti=N*B*abs(det(J))*trong_so(j);
            unit_vt(t(1:4,i))=unit_vt(t(1:4,i)) + unit_vti; % tao vecto don vi cho tich phan gauss
            
        end
    end
    

    %%
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

A_fix(I_active,I_active)=mtA(I_active,I_active);
mt_dA_dhp_fix(I_active,I_active)=mt_dA_dhp(I_active,I_active);
mt_dA_dalpha_r_fix(I_active,I_active)=mt_dA_dalpha_r(I_active,I_active);
mt_dA_dalpha_theta_fix(I_active,I_active)=mt_dA_dalpha_theta(I_active,I_active);

mt_dA_dhp_fix(I_active,I_boundry)=0;mt_dA_dhp_fix(I_boundry,I_active)=0;
mt_dA_dalpha_r_fix(I_active,I_boundry)=0;mt_dA_dalpha_r_fix(I_boundry,I_active)=0;
mt_dA_dalpha_theta_fix(I_active,I_boundry)=0;mt_dA_dalpha_theta_fix(I_boundry,I_active)=0;

b(I_boundry)=0;%still not fix but straightly to b
dB_dhp(I_boundry)=0;
dB_dalpha_r(I_boundry)=0;
dB_dalpha_theta(I_boundry)=0;

% bh(I_boundry)=0;
vtp=(A_fix\b)'/hp^2;% tinh vecto ap suat

Fz=unit_vt'*vtp';

M_X=vt_X'*vtp';
M_Y=vt_Y'*vtp';

Xc_tinh=M_X/Fz;
Yc_tinh=M_Y/Fz;

% dp_dhp=(A_fix\(dB_dhp-mt_dA_dhp_fix*vtp'));
dp_dalpha_r=(A_fix\(dB_dalpha_r-mt_dA_dalpha_r_fix*vtp'));
dp_dalpha_theta=(A_fix\(dB_dalpha_theta-mt_dA_dalpha_theta_fix*vtp'));

% dFz_dhp=unit_vt'*dp_dhp;
% dX_dhp=1/Fz.*(vt_X'*dp_dhp-dFz_dhp*XcF);
% dY_dhp=1/Fz.*(vt_Y'*dp_dhp-dFz_dhp*YcF);

dFz_dapha_r=unit_vt'*dp_dalpha_r;
dX_dapha_r=1/Fz.*(vt_X'*dp_dalpha_r-dFz_dapha_r*Xc_tinh);
dY_dapha_r=1/Fz.*(vt_Y'*dp_dalpha_r-dFz_dapha_r*Yc_tinh);

dFz_dapha_theta=unit_vt'*dp_dalpha_theta;
dX_dapha_theta=1/Fz.*(vt_X'*dp_dalpha_theta-dFz_dapha_theta*Xc_tinh);
dY_dapha_theta=1/Fz.*(vt_Y'*dp_dalpha_theta-dFz_dapha_theta*Yc_tinh);

% D=[dFz_dhp dFz_dapha_r dFz_dapha_theta;dX_dhp dX_dapha_r dX_dapha_theta;dY_dhp dY_dapha_r dY_dapha_theta];
% new_variables=[hp;alpha_r;alpha_theta]-D\[Fz-;XcF-X_p;YcF-Y_p]
% D=[dFz_dhp dFz_dapha_theta ;dX_dhp dX_dapha_theta ;dY_dhp dY_dapha_theta ];
% new_variables=[hp;alpha_theta]-D\[Fz-Wz;XcF-X_p;YcF-Y_p]
% D=[dFz_dhp dFz_dapha_r ;dX_dhp dX_dapha_r ;dY_dhp dY_dapha_r ];
% new_variables=[hp;alpha_r]-D\[Fz-Wz;XcF-X_p;YcF-Y_p]
% hp_new=new_variables(1);

D=[dX_dapha_r dX_dapha_theta;dY_dapha_r dY_dapha_theta];
new_variables=[alpha_r;alpha_theta]-D\[Xc_tinh-X_p;Yc_tinh-Y_p]
alpha_r_new=new_variables(1);
alpha_theta_new=new_variables(2);
(Xc_tinh-X_p)/Xc_tinh
(Yc_tinh-Y_p)/Yc_tinh

% if (norm(hp_new-hp)/norm(hp_new)<=1e-5&&...
% if (         norm(alpha_r_new-alpha_r)/norm(alpha_r_new)<=1e-5&&...
%         norm(alpha_theta_new-alpha_theta)/norm(alpha_theta_new)<=1e-5&&...
%         norm((Xc_tinh-X_p)/norm(Xc_tinh))<=1e-5&&...
%         norm((Yc_tinh-Y_p)/norm(Yc_tinh))<=1e-5)
%     disp('done')
% break
% end
% thu=thu+1;
%     end
mtp(:,:)=vtp(luoi(:,:));

mtp_bar(:,:)=vtp(luoi(:,:))*hp^2/(r1^2*omega*nuy0);
max(vtp)
min(H)*hp
%% toa do va luoi phan tu



Coordinate3D_polar=toadophantu3D(snut3D,kchia_m,kchia_n,kchia_q,r1,r2,beta,Y3D,H3D_global,luoi3D);
Coordinate3D_polar(3,:)=Coordinate3D_polar(3,:);
Coordinate3D_polar(2,:)=Coordinate3D_polar(2,:);
% plot3(Coordinate3D_polar(1,:),Coordinate3D_polar(2,:),Coordinate3D_polar(3,:),'-')
% plot(1:size(Coordinate3D_polar(2,:)),Coordinate3D_polar(2,:),'.')
% Coordinate3D_nature=zeros(3,snut3D);
% Coordinate3D_nature(1,:)=Coordinate3D_polar(2,:).*cos(Coordinate3D_polar(1,:));
% Coordinate3D_nature(3,:)=Coordinate3D_polar(3,:).*sin(Coordinate3D_polar(1,:));
% Coordinate3D_nature(2,:)=Coordinate3D_polar(2,:);
%% tinh dp/dtheta dp/dr (can kiem tra)
dp_dtheta_bar=zeros(snut,1);
dp_dr_bar=zeros(snut,1);
for i=1:kchia_m+1
    for j=2:kchia_n+1
        dp_dtheta_bar(luoi(i,j))=(mtp_bar(i,j)-mtp_bar(i,j-1))/(1/kchia_n);
    end
    dp_dtheta_bar(luoi(i,1))=interp1(2:kchia_n+1,dp_dtheta_bar(luoi(i,2:kchia_n+1)),1,'linear','extrap');
end

%%
for j=1:kchia_n+1
    for i=kchia_m+1:-1:2
        dp_dr_bar(luoi(i,j))=(mtp_bar(i-1,j)-mtp_bar(i,j))/((r2-r1)/kchia_m/r1);
    end
    dp_dr_bar(luoi(1,j))=interp1(2:kchia_m+1,dp_dr_bar(luoi(2:kchia_m+1,j)),1,'linear','extrap');
end

%% Tinh tich phan 1D
index_vertical=10;% so phan tu lay tich phan 1D
nuy_global(1:snut3D,1:kchia_q+1)=1;
Nuy_global_node(1:snut3D)=1;
nuy_global_chotichphan1D(1:snut3D,1:index_vertical+1)=1;%%%%%%%change later
unit_vt_1D_chotichphan1D=zeros(snut3D,index_vertical+1);
unit_vt_y_1D_chotichphan1D=zeros(snut3D,index_vertical+1);
Y3D_bar=Y3D(:)./H3D_global(:);
H3D_global_bar=H3D_global;
PHUONGCHIEUCAO=zeros(snut3D,kchia_q);% phuong chieu cao tai moi diem lay tich phan 1D
for i=1:snut3D
    if i>snut
    PHUONGCHIEUCAO(i,1:index_vertical+1)=Y3D_bar(i)/index_vertical*(0:index_vertical);% phuong chieu cao theo do cao 3D
    nuy_local=interp1(H3D_global_bar(i)/kchia_q*(0:kchia_q),nuy_global(i,1:kchia_q+1),PHUONGCHIEUCAO(i,1:index_vertical+1),'linear','extrap'); %noi suy 
    for k=1:index_vertical %chi so chieu cao
        for j=1:2
            [xi,trong_so]=tpso_Gauss_1d(2,j);
            [N1,N2]=shapefunc2D(xi);
            N=[N1 N2];
            L=Y3D_bar(i)/index_vertical;%det J
            bien_y_i=N*[PHUONGCHIEUCAO(i,k);PHUONGCHIEUCAO(i,k+1)];
            
            nuy_1D_i=N*[nuy_local(k);nuy_local(k+1)];
            
            unit_vt_1D_i=N*L/2*trong_so(j);
            unit_vt_1D_chotichphan1D(i,k)=unit_vt_1D_chotichphan1D(i,k)+unit_vt_1D_i(1);
            unit_vt_1D_chotichphan1D(i,k+1)=unit_vt_1D_chotichphan1D(i,k+1)+unit_vt_1D_i(2);
            
            unit_vt_y_1D_i=N*bien_y_i*L/2*trong_so(j);
            unit_vt_y_1D_chotichphan1D(i,k)=unit_vt_y_1D_chotichphan1D(i,k)+unit_vt_y_1D_i(1);
            unit_vt_y_1D_chotichphan1D(i,k+1)=unit_vt_y_1D_chotichphan1D(i,k+1)+unit_vt_y_1D_i(2);
            
        end
    end
    end
end
%% Tinh cac F
F0=zeros(snut3D,1);
F1=zeros(snut3D,1);

for i =1:snut3D
        F0(i)=1./nuy_global_chotichphan1D(i,:)*unit_vt_1D_chotichphan1D(i,:)';
        F1(i)=1./nuy_global_chotichphan1D(i,:)*unit_vt_y_1D_chotichphan1D(i,:)';
end



%% chuyen doi medium to max
chiso_max=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
chiso_min=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
for k=1:kchia_q+1
    for j=1:kchia_n+1
        for i=1:kchia_m+1
            chiso_max(luoi3D(i,j,k))=luoi3D(i,j,kchia_q+1);
            chiso_min(luoi3D(i,j,k) )=luoi3D(i,j,1);
            
        end
    end
end

%% tinh van toc u w
u1_bar=zeros(snut3D,1);
for i=1:kchia_n+1
    for j=1:kchia_q+1
%         u1_bar(luoi3D(kchia_m+1:-1:1,i,j))=(r1+(0:kchia_m)/kchia_m*(r2-r1))/r1;
u1_bar(luoi3D(kchia_m+1:-1:1,i,j))=1;
    end
end
A1=zeros(snut3D,1);
A2=zeros(snut3D,1);
B1=zeros(snut3D,1);
B2=zeros(snut3D,1);
for i=1:snut3D
    r=Coordinate3D_polar(3,i);
    A1(i)=(H3D_global(i))^2*dp_dr_bar(chiso_min(i))*(F1(i)-F1(chiso_max(i))/F0(chiso_max(i))*F0(i));
    A2(i)=(H3D_global(i))^2/(beta*r^2)*dp_dtheta_bar(chiso_min(i))*(F1(i)-F1(chiso_max(i))/F0(chiso_max(i))*F0(i))-1/F0(chiso_max(i))*F0(i)+1;
    B1(i)=u1_bar(i)/F0(chiso_max(i))*F0(i);
    B2(i)=0;
    
end
vel_u_bar=zeros(snut3D,1);
vel_w_bar=zeros(snut3D,1);
vel_v_bar=zeros(snut3D,1);
for i=1:snut3D
    vel_w_bar(i)=A1(i);
%     vel_u_bar(i)=vel(1);
    vel_u_bar(i)=A2(i);
%     vel_v(i)=0;%%%%%%%%%%%%%%%%change later
end


%% tinh du/dx dw/dr (can kiem tra)
% du_dtheta_bar=zeros(snut3D,1);
% dw_dr_bar=zeros(snut3D,1);
du_dy_bar=zeros(snut3D,1);
dw_dy_bar=zeros(snut3D,1);
% for k=1:kchia_q+1
%     for i=1:kchia_m+1
%         for j=2:kchia_n+1
%             du_dtheta_bar(luoi3D(i,j,k))=(vel_u_bar(luoi3D(i,j,k))-vel_u_bar(luoi3D(i,j-1,k)))/beta*kchia_n;
%         end
%         du_dtheta_bar(luoi3D(i,1,k))=interp1(2:kchia_n+1,du_dtheta_bar(luoi3D(i,2:kchia_n+1,k)),1,'linear','extrap');
%     end
% end


% for k=1:kchia_q+1
%     for j=1:kchia_n+1
%         for i=kchia_m+1:-1:2
%             dw_dr_bar(luoi3D(i,j,k))=(vel_w_bar(luoi3D(i-1,j,k))-vel_w_bar(luoi3D(i,j,k)))/(r2-r1)*kchia_m*r1;
%         end
%         dw_dr_bar(luoi3D(1,j,k))=interp1(2:kchia_m+1,dw_dr_bar(luoi3D(2:kchia_m+1,j,k)),1,'linear','extrap');
%     end
% end

du_dtheta_bar=d_dtheta_bar(vel_u_bar);
dw_dr_bar=d_dr_bar(vel_w_bar);
%d/dy
for i=1:kchia_m+1
    for j=1:kchia_n+1
        for k=2:kchia_q+1
            du_dy_bar(luoi3D(i,j,k))=(vel_u_bar(luoi3D(i,j,k))-vel_u_bar(luoi3D(i,j,k-1)))*kchia_q;
            dw_dy_bar(luoi3D(i,j,k))=(vel_w_bar(luoi3D(i,j,k))-vel_w_bar(luoi3D(i,j,k-1)))*kchia_q;
        end
        du_dy_bar(luoi3D(i,j,1))=interp1(2:kchia_q+1,du_dy_bar(luoi3D(i,j,2:kchia_q+1)),1,'linear','extrap');
        dw_dy_bar(luoi3D(i,j,1))=interp1(2:kchia_q+1,dw_dy_bar(luoi3D(i,j,2:kchia_q+1)),1,'linear','extrap');
    end
end
%% tinh vel_v
vel_v_unit=zeros(snut3D,index_vertical+1);
vel_v_tp1=zeros(snut3D,index_vertical+1);
vel_v_tp2=zeros(snut3D,index_vertical+1);
du_dtheta_global_bar=convert3Dfor1D(kchia_m,kchia_n,kchia_q,luoi3D,du_dtheta_bar);
dw_dr_global_bar=convert3Dfor1D(kchia_m,kchia_n,kchia_q,luoi3D,dw_dr_bar);
h_bar=Y3D(:)./H3D_global(:);
for i=1:snut3D
    if i>snut
    PHUONGCHIEUCAO(i,1:index_vertical+1)=h_bar(i)/index_vertical*(0:index_vertical);% phuong chieu cao theo do cao 3D
%     nuy_local=interp1(H3D_global_bar(i)/kchia_q*(0:kchia_q),nuy_global(i,1:kchia_q+1),PHUONGCHIEUCAO(i,1:index_vertical+1),'linear','extrap'); %noi suy 
%     du_dtheta_local=interp1(H3D_global_bar(i)/kchia_q*(0:kchia_q),du_dtheta_global_bar(i,1:kchia_q+1),PHUONGCHIEUCAO(i,1:index_vertical+1),'linear','extrap'); %noi suy 
%     dw_dr_local=interp1(H3D_global_bar(i)/kchia_q*(0:kchia_q),dw_dr_global_bar(i,1:kchia_q+1),PHUONGCHIEUCAO(i,1:index_vertical+1),'linear','extrap'); %noi suy 
    for k=1:index_vertical %chi so chieu cao
        for j=1:2
            [xi,trong_so]=tpso_Gauss_1d(2,j);
            [N1,N2]=shapefunc2D(xi);
            N=[N1 N2];
            L=h_bar(i)/index_vertical;%det J
            bien_y_i=N*[PHUONGCHIEUCAO(i,k);PHUONGCHIEUCAO(i,k+1)];
            
%             nuy_1D_i=N*[nuy_local(k);nuy_local(k+1)];
            
%             N*(du_dtheta(k:k+1)+dw_dr(k:k+1))'
%             (N*(du_dtheta_global(i,k:k+1)+dw_dr_global(i,k:k+1))');
            vel_v_i=-N*H3D_global_bar(i)*(N*(du_dtheta_global_bar(i,k:k+1)+dw_dr_global_bar(i,k:k+1))')*L/2*trong_so(j);
            vel_v_unit(i,k)=vel_v_unit(i,k)+vel_v_i(1);
            vel_v_unit(i,k+1)=vel_v_unit(i,k+1)+vel_v_i(2);
            
            vel_v_tp1_i=N*Coordinate3D_polar(3,i)*H3D_global(i)*vel_w_bar(i)*L/2*trong_so(j);
            vel_v_tp2_i=N*H3D_global(i)*vel_u_bar(i)*L/2*trong_so(j);
            
            vel_v_tp1(i,k)=vel_v_tp1(i,k)+vel_v_tp1_i(1);
            vel_v_tp2(i,k+1)=vel_v_tp2(i,k+1)+vel_v_tp2_i(2);
%             unit_vt_y_1D_i=N*bien_y_i*L/2*W(j);
%             unit_vt_y_1D(i,k)=unit_vt_y_1D(i,k)+unit_vt_y_1D_i(1);
%             unit_vt_y_1D(i,k+1)=unit_vt_y_1D(i,k+1)+unit_vt_y_1D_i(2);
        end
    end
    end
    vel_v_bar(i)=sum(vel_v_unit(i,:));
    vel_v_tp1_bar(i)=sum(vel_v_tp1(i,:));
    vel_v_tp2_bar(i)=sum(vel_v_tp2(i,:));
end
%%
d_tp1_r=d_dr_bar(vel_v_tp1_bar);
d_tp2_theta=d_dtheta_bar(vel_v_tp2_bar);
dh_dtheta_bar=d_dtheta_bar(H3D_global);
dh_dr_bar=d_dr_bar(H3D_global);
for i=1:snut3D
   vel_v_bar(i)=-1/Coordinate3D_polar(3,i)*d_tp1_r(i)-1/beta*d_tp2_theta(i)+Coordinate3D_polar(2,i)*dh_dr_bar(i)*vel_w_bar(i)+Coordinate3D_polar(2,i)/beta*dh_dtheta_bar(i)*vel_u_bar(i);
end
%% tinh nhiet


%%
node_index3D=chuyennutvaophantu3D(spt3D,luoi_phantu3D,luoi3D);
%% tich phan nhiet do
K=zeros(snut_all,snut_all);
K1_1=zeros(snut_all,snut_all);
K1_2=zeros(snut_all,snut_all);
f=zeros(snut_all,1);
%he so
f1=k_conduct/(delta*Cp*omega*hp^2);
f2=nuy0*r1^2*omega/(delta*Cp*T0*hp^2);
for i=1:spt3D
    
    xP=Coordinate3D_polar(1,node_index3D(1:8,i));
    yP=Coordinate3D_polar(3,node_index3D(1:8,i));
    zP=Coordinate3D_polar(2,node_index3D(1:8,i));
    
%     xP_bar=Coordinate3D_polar(1,node_index3D(1:8,i))/r1;
%     yP_bar=Coordinate3D_polar(3,node_index3D(1:8,i))./H3D_global(node_index3D(1:8,i));
%     zP_bar=Coordinate3D_polar(2,node_index3D(1:8,i))/r1;
    
    for j=1:8
        [R,S,T,trong_so]=tpso_Gauss_3d(2,j);
        [N1,N2,N3,N4,N5,N6,N7,N8]=hamdang3D(R,S,T);
        N=[N1,N2,N3,N4,N5,N6,N7,N8]';
        dNrst=dNrst_3D(R,S,T);
        J3D=Jctugiac3D(xP,yP,zP,R,S,T);
        Jn=inv(J3D);
        dNx(1:8)=Jn(1,1)*dNrst(1:8,1)+Jn(1,2)*dNrst(1:8,2)+Jn(1,3)*dNrst(1:8,3);
        dNy(1:8)=Jn(2,1)*dNrst(1:8,1)+Jn(2,2)*dNrst(1:8,2)+Jn(2,3)*dNrst(1:8,3);
        dNz(1:8)=Jn(3,1)*dNrst(1:8,1)+Jn(3,2)*dNrst(1:8,2)+Jn(3,3)*dNrst(1:8,3);
        
        A=N'*xP';%phuong theta
        B=N'*yP';%phuong chieu cao
        C=N'*zP';%phuong ban kinh
        U_=N'*vel_u_bar(node_index3D(1:8,i));
        V_=N'*vel_v_bar(node_index3D(1:8,i));
        W_=N'*vel_w_bar(node_index3D(1:8,i));
        HH=N'*(H3D_global_bar(node_index3D(1:8,i)));
        DHX=dNx*(H3D_global_bar(node_index3D(1:8,i)));
        DHZ=dNz*(H3D_global_bar(node_index3D(1:8,i)));
        Nuy=N'*(Nuy_global_node(node_index3D(1:8,i)))';
        
        K_i=N*(1/beta*U_*dNx+1/HH*(V_-beta*U_*B*DHX-W_*B*DHZ)*dNy+W_*dNz+W_/C*N')*abs(det(J3D))*C*trong_so(j);
        K(node_index3D(1:8,i),node_index3D(1:8,i))=K(node_index3D(1:8,i),node_index3D(1:8,i))+K_i;
        K1_i=f1/HH^2*(dNy'*dNy)*abs(det(J3D))*C*trong_so(j);% ++++++++
        K1_1(node_index3D(1:8,i),node_index3D(1:8,i))=K1_1(node_index3D(1:8,i),node_index3D(1:8,i))+K1_i;
        
        Du_dy=N'*(du_dy_bar(node_index3D(1:8,i)));
        Dw_dy=N'*(dw_dy_bar(node_index3D(1:8,i)));
        
        f_i=N*f2*Nuy/HH^2*(C^2*Du_dy^2+Dw_dy^2)*abs(det(J3D))*C*trong_so(j);
        f(node_index3D(1:8,i))=f(node_index3D(1:8,i)) + f_i;
    end

end
% plot3(xP,yP,zP,'-')
% tinh K1_2

% for i=1:spt2D
%     theta=p(1,t(1:4,i));
%     r=p(2,t(1:4,i));
%     xP=pp(1,t(1:4,i));
%     yP=pp(2,t(1:4,i));
%     for j=1:4
%         [xi,eta,trong_so]=tpso_Gauss_2d(2,j);
%         [N1,N2,N3,N4]=hamdang(xi,eta);
%         N=[N1 ;N2 ;N3 ;N4];
%         dNij=[ -1/4*(1-eta) -1/4*(1-xi);...
%             1/4*(1-eta) -1/4*(1+xi);...
%             1/4*(1+eta)  1/4*(1+xi);...
%             -1/4*(1+eta)  1/4*(1-xi)];
%         
%         J=Jctugiac(theta,r,eta,xi);
%         Jn=inv(J);
%         
%         dNx(1:4)=Jn(1,1)*dNij(1:4,1)+Jn(1,2)*dNij(1:4,2);
%         dNz(1:4)=Jn(2,1)*dNij(1:4,1)+Jn(2,2)*dNij(1:4,2);
%         
%         A=N'*theta';%bien
%         B=N'*r';%bien
%         C=N'*xP';%toa do de cac
%         D=N'*yP';%toa do de cac\
%         
%         HH=N'*(H(t(1:4,i)))';
%         HH2=HH^2;
%         
% %         K1_2_i=N*(u_bar
%     end
% end

%% tinh K nhiet pad



K_pad1=zeros(snut_all);
node_index3D_pad_all=chuyennutvaophantu3D(spt3D*2,luoi_phantu3D_all,luoi_all);
khoang_i_pad_node=min(luoi_pad3D,[],'all'):max(luoi_pad3D,[],'all');
khoang_i_pad_phantu=(spt3D+1):1:(2*spt3D);
Coordinate3D_polar_pad_bar=toadophantu3D_pad(snut3D,kchia_m,kchia_n,kchia_q,r1,r2,beta,luoi_pad3D,khoang_i_pad_node,t_pad);
% plot3(Coordinate3D_polar_pad_bar(1,:),Coordinate3D_polar_pad_bar(2,:),Coordinate3D_polar_pad_bar(3,:),'.')
for i=khoang_i_pad_phantu
    
    xP=Coordinate3D_polar_pad_bar(1,node_index3D_pad_all(1:8,i));
    yP=Coordinate3D_polar_pad_bar(2,node_index3D_pad_all(1:8,i));
    zP=Coordinate3D_polar_pad_bar(3,node_index3D_pad_all(1:8,i));
    
%     xP_bar=Coordinate3D_polar(1,node_index3D(1:8,i))/r1;
%     yP_bar=Coordinate3D_polar(3,node_index3D(1:8,i))./H3D_global(node_index3D(1:8,i));
%     zP_bar=Coordinate3D_polar(2,node_index3D(1:8,i))/r1;
    
    for j=1:8
        [R,S,T,trong_so]=tpso_Gauss_3d(2,j);
        [N1,N2,N3,N4,N5,N6,N7,N8]=hamdang3D(R,S,T);
        N=[N1,N2,N3,N4,N5,N6,N7,N8]';
        dNrst=dNrst_3D(R,S,T);
        J3D=Jctugiac3D(xP,yP,zP,R,S,T);
        Jn=inv(J3D);
        dNx(1:8)=Jn(1,1)*dNrst(1:8,1)+Jn(1,2)*dNrst(1:8,2)+Jn(1,3)*dNrst(1:8,3);
        dNy(1:8)=Jn(2,1)*dNrst(1:8,1)+Jn(2,2)*dNrst(1:8,2)+Jn(2,3)*dNrst(1:8,3);
        dNz(1:8)=Jn(3,1)*dNrst(1:8,1)+Jn(3,2)*dNrst(1:8,2)+Jn(3,3)*dNrst(1:8,3);
        
        A=N'*xP';%phuong theta
        B=N'*yP';%phuong chieu cao
        C=N'*zP';%phuong ban kinh
        
                
        Kpad1_i=(dNx'*dNx/B^2+dNy'*dNy*B+dNz'*dNz*B-N*dNy)*C*abs(det(J3D))*trong_so(j);
        K_pad1(node_index3D_pad_all(1:8,i),node_index3D_pad_all(1:8,i))=K_pad1(node_index3D_pad_all(1:8,i),node_index3D_pad_all(1:8,i))+Kpad1_i;
        
        
    end

end
%% tinh cac tich phan bien
% z1
luoi_bien2D_z1=squeeze(luoi_pad3D(:,:,kchia_q+1));
luoi_phantu_bien2D_z1=luoiphantu2D(kchia_m,kchia_n);
node_index_bien2D_z1=chuyennutvaophantu(kchia_m*kchia_n,luoi_phantu_bien2D_z1,kchia_m,kchia_n,luoi_bien2D_z1);
toado_z1=Coordinate3D_polar_pad_bar([1,2],luoi_all(:));
[mtA_z1,f_z1]=tichphanbien2D(node_index_bien2D_z1, kchia_m*kchia_n,toado_z1, T0, snut_all,3);
% % mat tiep xuc
% luoi_bien2D_z2=squeeze(luoi_pad3D(:,:,1));
% luoi_phantu_bien2D_z2=luoiphantu2D(kchia_m,kchia_n);
% node_index_bien2D_z2=chuyennutvaophantu(kchia_m*kchia_n,luoi_phantu_bien2D_z2,kchia_m,kchia_n,luoi_bien2D_z2);
%r1
luoi_bien2D_r1=squeeze(luoi_pad3D(kchia_m+1,:,:));
luoi_phantu_bien2D_r1=luoiphantu2D(kchia_n,kchia_q);
node_index_bien2D_r1=chuyennutvaophantu(kchia_n*kchia_q,luoi_phantu_bien2D_r1,kchia_n,kchia_q,luoi_bien2D_r1);
toado_r1=Coordinate3D_polar_pad_bar([1,3],luoi_all(:));
[mtA_r1,f_r1]=tichphanbien2D(node_index_bien2D_r1, kchia_m*kchia_q,toado_r1, T0, snut_all,1);
%r2
luoi_bien2D_r2=squeeze(luoi_pad3D(1,:,:));
luoi_phantu_bien2D_r2=luoiphantu2D(kchia_n,kchia_q);
node_index_bien2D_r2=chuyennutvaophantu(kchia_n*kchia_q,luoi_phantu_bien2D_r2,kchia_n,kchia_q,luoi_bien2D_r2);
toado_r2=Coordinate3D_polar_pad_bar([1,3],luoi_all(:));
[mtA_r2,f_r2]=tichphanbien2D(node_index_bien2D_r1, kchia_m*kchia_q,toado_r2, T0, snut_all,2);
% fi1
luoi_bien2D_fi1=squeeze(luoi_pad3D(:,kchia_n+1,:));
luoi_phantu_bien2D_fi1=luoiphantu2D(kchia_m,kchia_q);
node_index_bien2D_fi1=chuyennutvaophantu(kchia_m*kchia_q,luoi_phantu_bien2D_fi1,kchia_m,kchia_q,luoi_bien2D_fi1);
toado_fi1=Coordinate3D_polar_pad_bar([2,3],luoi_all(:));
[mtA_fi1,f_fi1]=tichphanbien2D(node_index_bien2D_fi1, kchia_m*kchia_q,toado_fi1, T0, snut_all,4);
% fi2
luoi_bien2D_fi2=squeeze(luoi_pad3D(:,1,:));
luoi_phantu_bien2D_fi2=luoiphantu2D(kchia_m,kchia_q);
node_index_bien2D_fi2=chuyennutvaophantu(kchia_m*kchia_q,luoi_phantu_bien2D_fi2,kchia_m,kchia_q,luoi_bien2D_fi2);
toado_fi2=Coordinate3D_polar_pad_bar([2,3],luoi_all(:));
[mtA_fi2,f_fi2]=tichphanbien2D(node_index_bien2D_fi2, kchia_m*kchia_q,toado_fi2, T0, snut_all,5);

%%  Ghep matran do cung
K_all=K+K1_1+K_pad1+mtA_z1+mtA_r1+mtA_r2+mtA_fi1+mtA_fi2;
boundaryinlet_oil1=squeeze(luoi3D(:,1,:));
boundaryinlet_oil2=squeeze(luoi3D(1,:,:));
boundaryinlet_oil3=squeeze(luoi3D(:,:,1));
boundaryinlet_oil4=squeeze(luoi3D(:,kchia_n+1,:));
boundaryinlet_oil5=squeeze(luoi3D(kchia_m+1,:,:));
boundaryinlet_oil6=squeeze(luoi3D(:,:,kchia_q+1));
% boundaryinlet_oil=[boundaryinlet_oil1(:),boundaryinlet_oil2(:),boundaryinlet_oil3(:),boundaryinlet_oil4(:),boundaryinlet_oil5(:)];
f_all=f-f_z1-f_r1+f_r2-f_fi1+f_fi2;
K_all(boundaryinlet_oil1,:)=0;K_all(:,boundaryinlet_oil1)=0;
% (K_all(1:112,1:112))
K_oil=K+K1_1;
f_oil=f;
K_oil(boundaryinlet_oil1,:)=0;
% K_oil(:,boundaryinlet_oil)=0;

for i=1:length(boundaryinlet_oil1(:))
    j=boundaryinlet_oil1(i);
    K_all(j,j)=1;
    K_oil(j,j)=1;
end
% (K_all(1:112,1:112))
f_oil(boundaryinlet_oil1)=T0/T0;
%% Tinh nhiet
% K_part1=K_all(1:snut3D,1:snut3D);
% T=K_all\f_all;
T1=K_oil(1:snut3D,1:snut3D)\f_oil(1:snut3D);
T1(boundaryinlet_oil6)
%% update
%%
    % cd=linspace(0,1,kchia_m+1);
    % mtp(:,:)=vtp(luoi(:,:));
    % mtq(:,:)=vtq(luoi(:,:));
    % [X,Y]=meshgrid(cd,phi);
    % figure(1)
    % surf(X',Y',mtp)
    % figure(2)
    % surf(X',Y',mtq)
    
    
    
    % vtp_giua=vtp(luoi(kchia_m/2,:));
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
%% ve hinh
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
figure(3)
plot(0:kchia_m,mtp(kchia_m/2+1,1:kchia_n+1),0:kchia_m,mtp(kchia_m*0.8+1,1:kchia_n+1),0:kchia_m,mtp(kchia_m*0.2+1,1:kchia_n+1))
figure(4)
plot(kchia_n:-1:0,mtp(1:kchia_m+1,kchia_n*0.7+1),kchia_n:-1:0,mtp(1:kchia_m+1,kchia_n*0.6+1),kchia_n:-1:0,mtp(1:kchia_m+1,kchia_n*0.1+1))



%%
% figure(5)
% plot3(x,y,du_dtheta_bar,'.')
% 
% figure(6)
% plot3(x,y,dp_dr_bar,'.')