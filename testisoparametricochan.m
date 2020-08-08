clc
clear all

thu=1;
kchia_m=40;kchia_n=40;kchia_q=20;%so khoang chia
snut=(kchia_m+1)*(kchia_n+1);
snut3D=(kchia_m+1)*(kchia_n+1)*(kchia_q+1);
spt2D=kchia_m*kchia_n;
spt3D=kchia_m*kchia_n*kchia_q;

% beta=20/360*2*pi;
% r1=1.105;
% r2=1.690;
% nuy=0.004043;%do nhot dong luc hoc (N/m^2.s)
% omega=107/60*(2*pi);%toc do goc (rad/s)
% h2=0.000005;
% h1=0.000015;

% beta=55/360*2*pi;
% r1=0.02;
% r2=0.05;
% nuy=0.039;%do nhot dong luc hoc (N/m^2.s)
% omega=3000/60*(2*pi);%toc do goc (rad/s)

beta=50/360*2*pi;
r1=0.05715;
r2=0.1143;
nuy=0.039;%do nhot dong luc hoc (N/m^2.s)
omega=1500/60*(2*pi);%toc do goc (rad/s)

%% danh so nut
p_so=1;
luoi=zeros(kchia_m+1,kchia_n+1);
luoi3D=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
chiso_max=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
chiso_min=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
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
%% luoi 3D
p_so=1;
for k=1:kchia_q+1
    for j=1:kchia_n+1
        if mod(j,2)~=0
            for i=1:kchia_m+1
                luoi3D(i,j,k)=p_so;
                p_so=p_so+1;
            end
        end
        if mod(j,2)==0
            for i=1:kchia_m+1
                luoi3D(i,j,k)=p_so;
                p_so=p_so+1;
            end
        end
    end
end
%% chuyen doi medium to max
for k=1:kchia_q+1
    for j=1:kchia_n+1
        for i=1:kchia_m+1
            chiso_max(luoi3D(i,j,k))=luoi3D(i,j,kchia_q+1);
            chiso_min(luoi3D(i,j,k))=luoi3D(i,j,1);
            
        end 
    end
end

%% danh so phan tu
p_so2=1;
luoi_phantu=zeros(kchia_m,kchia_n);
luoi_phantu3D=zeros(kchia_m,kchia_n,kchia_q);
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
%3D
dem=1;
for k=1:kchia_q
    for j=1:kchia_n
        for i=1:kchia_m
            luoi_phantu3D(i,j,k)=dem;
            dem=dem+1;
        end
    end
end
%% toa do va luoi phan tu

p=zeros(2,snut);
for i=1:snut
    [row,col]=find(luoi==i);
    p(1,i)=(1-col)/kchia_n*L;
    p(2,i)=(1/2-(kchia_m-row+1)/kchia_m)*L/(2*R)*2;
end


t=zeros(4,spt2D);
for i=1:kchia_m*kchia_n
    [row,col]=find(luoi_phantu==i);
%     if (col~=n)
        t(1:4,i)=[luoi((row+1),col),luoi((row+1),col+1),luoi(row,col+1),luoi(row,col)];
%     else
%         t(1:4,i)=[luoi((row+1),col),luoi((row+1),1),luoi(row,1),luoi(row,col)];
%     end
end


%%
    Wz=52265;
%     while(1)
        if (thu==1)

        h2=0.000075;
        h1=h2*2.75;
        deltah=h1-h2;
        else

        h2=h2_new;
        h1=h2+deltah;
            
        end
        
        delta=(h1-h2)/h2;
        
    dem_ib=1;
    dem_ia=1;
    h=zeros(kchia_m+1,kchia_n+1);
%     phuongchieucao=zeros(kchia_n+1,kchia_n+1,kchia_n+1);%%%%%%%%%%%%%%%%%%
%     PHUONGCHIEUCAO=zeros((kchia_n+1)*(kchia_n+1),kchia_n+1);%%%%%%%%%%%%%%%%%%
    for i=1:kchia_m+1
        for j=1:kchia_n+1
            h(i,j)=1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/beta);
            %h(i,j)=(1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/sqrt(pp(1,luoi(i,j))^2 + pp(2,luoi(i,j))^2)));
            
%             phuongchieucao(i,j,1:kchia_q+1)=(0:kchia_q)*h(i,j);%%%%%%%%%%%%%%%%
%             PHUONGCHIEUCAO(luoi(i,j),1:kchia_q+1)=(0:kchia_q)*h(i,j);
            if (i==1||j==1||i==kchia_m+1||j==kchia_n+1)
            I_boundry(dem_ib)=luoi(i,j);
            dem_ib=dem_ib+1;
            else
            I_active(dem_ia)=luoi(i,j);
            dem_ia=dem_ia+1;
            end
        end
    end
%%    
    H(luoi(:,:))=h(:,:);% chuyen array 2d h thanh H 1d
    H3D=zeros((kchia_m+1)*(kchia_n+1)*(kchia_q+1),1);
    H3D_global=zeros((kchia_m+1)*(kchia_n+1)*(kchia_q+1),1);
    luoi_max=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
    for i=1:kchia_m+1
        for j=1:kchia_n+1
            for k=2:kchia_q+1
                H3D(luoi3D(i,j,k))=h(i,j)/(kchia_q+2-k);
                H3D_global(luoi3D(i,j,k))=h(i,j);
                
            end
        end
    end
    
%%  
    % tich phan so
    mtA=zeros(snut,snut);
    mtA_d=zeros(snut,snut);
    mtAh=zeros(snut,snut);
    mtAh_d=zeros(snut,snut);
    dmtA_dphi=zeros(snut,snut);
    dmtA_dr=zeros(snut,snut);
    dA_dphi=zeros(snut,snut);
    dA_dr=zeros(snut,snut);
    b=zeros(snut,1);
    bh=zeros(snut,1);
    db_dphi=zeros(snut,1);
    db_dr=zeros(snut,1);
    unit_vt=zeros(snut,1);
    vt_X=zeros(snut,1);
    vt_Y=zeros(snut,1);
    %tich phan so Gauss
    

    for i=1:spt2D
            xP=p(1,t(1:4,i));
            yP=p(2,t(1:4,i));
            
        for j=1:9
            [xi,eta,W]=tpso_Gauss_2d(3,j);
            [N1,N2,N3,N4]=hamdang(xi,eta);
            N=[N1 ;N2 ;N3 ;N4];
            dNij=[ -1/4*(1-eta) -1/4*(1-xi);...
                1/4*(1-eta) -1/4*(1+xi);...
                1/4*(1+eta)  1/4*(1+xi);...
                -1/4*(1+eta)  1/4*(1-xi)];
            
            J=Jctugiac(xP,yP,eta,xi);
            Jn=inv(J);
            
            dNx(1:4)=Jn(1,1)*dNij(1:4,1)+Jn(1,2)*dNij(1:4,2);
            dNz(1:4)=Jn(2,1)*dNij(1:4,1)+Jn(2,2)*dNij(1:4,2);
            
            A=N'*xP';
            B=N'*yP';

%             size(r)
%             HH3=(N'*H(t(1:4,i))')^3;
%             HH2=(N'*H(t(1:4,i))')^2;
%             HH=(N'*H(t(1:4,i))');
            HH=1 + delta*(1 - acot(B/A)/beta);
            HH2=(1 + delta*(1 - acot(B/A)/beta))^2;
            HH3=(1 + delta*(1 - acot(B/A)/beta))^3;
            %tinh cac ma tran do cung
            Ai=HH3*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
            Ah=3*HH2*-(beta - acot(B/A))*h1/(h2^2*beta)*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
            dA_dphi_i=3*HH2*-(h1 - h2)*B/(A^2*(B^2/A^2 + 1)*beta*h2)*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
            dA_dr_i=3*HH2*(h1 - h2)/(A*(B^2/A^2 + 1)*beta*h2)*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
            
            mtA(t(1:4,i),t(1:4,i))=mtA(t(1:4,i),t(1:4,i))+Ai;
            mtAh(t(1:4,i),t(1:4,i))=mtAh(t(1:4,i),t(1:4,i))+Ah;
            dA_dphi(t(1:4,i),t(1:4,i))=dA_dphi(t(1:4,i),t(1:4,i))+dA_dphi_i;
            dA_dr(t(1:4,i),t(1:4,i))=dA_dr(t(1:4,i),t(1:4,i))+dA_dr_i;
            %tinh vecto tai
%             A=xP;
%             B=yP;
%             dH=h2*delta*(B./(A.^2.*(B.^2./A.^2 + 1).*sqrt(A.^2 + B.^2)) + acot(B./A).*A./(A.^2 + B.^2).^(3/2));
            
            
            r=sqrt(A^2+B^2);
            

%             dH=(h2-h1)/h2*B/(A^2*(B^2/A^2 + 1)*beta);
            
            %dH=delta*(-B/(A^2*(B^2/A^2 + 1)*sqrt(A^2 + B^2)) + acot(B/A)*A/(A^2 + B^2)^(3/2));
            
            
            %Bj=6*N*(N'*dH')*(N'*sqrt(A.^2+B.^2)')*abs(det(J))*W(j);
%            Bj=-6*nuy*omega/h2^2*N*dH*r*abs(det(J))*W(j);
%            Bj=-6*nuy*omega/h2^2*N*-(h1 - h2)*B^2/((A^2 + B^2)*h2*beta)*abs(det(J))*W(j);
            %Bj=-6*nuy*omega/h2^2*N*r*-(h1 - h2)*B/(h2*(A^2 + B^2)*beta)*abs(det(J))*W(j);
            Bj=-6*nuy*omega/h2^2*N*-(h1 - h2)*B/(sqrt(A^2 + B^2)*h2*beta)*abs(det(J))*W(j)*sin(C)...
            +6*nuy*omega/h2^2*N*(h1 - h2)*A/(sqrt(A^2 + B^2)*beta*h2)*abs(det(J))*W(j)*cos(C);
            
            %bhi=-6*nuy*omega*N*(A/sqrt(A^2 + B^2))*abs(det(J))*W(j)/h2^2;
            bhi=-6*nuy*omega*N*B*(3*h1 - 2*h2)/(sqrt(A^2 + B^2)*h2^4*beta)*abs(det(J))*W(j);
            db_dphi_i=-6*nuy*omega*N*(h1 - h2)*B*A/((A^2 + B^2)^(3/2)*h2^3*beta)*abs(det(J))*W(j);
            db_dr_i=-6*nuy*omega*N*-(h1 - h2)*A^2/((A^2 + B^2)^(3/2)*h2^3*beta)*abs(det(J))*W(j);
            
            b(t(1:4,i))=b(t(1:4,i))+Bj;
            bh(t(1:4,i))=bh(t(1:4,i))+bhi;
            db_dphi(t(1:4,i))=db_dphi(t(1:4,i))+db_dphi_i;
            db_dr(t(1:4,i))=db_dr(t(1:4,i))+db_dr_i;
            
            unit_vti=N*abs(det(J))*W(j);
            unit_vt(t(1:4,i))=unit_vt(t(1:4,i)) + unit_vti;
            
            X=N*A*abs(det(J))*W(j);
            Y=N*B*abs(det(J))*W(j);
            
            vt_X(t(1:4,i))=vt_X(t(1:4,i)) + X;%vecto x cho chuyen doi nodal forces
            vt_Y(t(1:4,i))=vt_Y(t(1:4,i)) + Y;%vecto y cho chuyen doi nodal forces
            
            
        end
    end
    
for i=1:snut
    if (sum(I_active==i)==0) % i khong thuoc I active
        mtA_d(i,i)=1;
    end
end

mtA_d(I_active,I_active)=mtA(I_active,I_active);
mtAh_d(I_active,I_active)=mtAh(I_active,I_active);
dmtA_dphi(I_active,I_active)=dA_dphi(I_active,I_active);
dmtA_dr(I_active,I_active)=dA_dr(I_active,I_active);

b(I_boundry)=0;
bh(I_boundry)=0;
db_dphi(I_boundry)=0;
db_dr(I_boundry)=0;

vtp=(mtA_d\b)';
Fz=unit_vt'*vtp';

M_X=vt_X'*vtp';
M_Y=vt_Y'*vtp';

Xc_tinh=M_X/Fz;
Yc_tinh=M_Y/Fz;
%%
dp_d=(mtA_d\(bh-mtAh_d*vtp'));
dp_dphi=(mtA_d\(db_dphi-dmtA_dphi*vtp'));
dp_dr=(mtA_d\(db_dr-dmtA_dr*vtp'));
D=unit_vt'*dp_d;
h2_new=h2-D\(Fz-Wz);
% if ((norm(h2_new-h2)/norm(h2_new))<=1e-10&&(norm(Fz-Wz)/norm(Wz))<1e-7)
%     disp('done')
% break
% end
% thu=thu+1;
%     end

%% Tinh tich phan 1D
index_vertical=6;% so phan tu lay tich phan 1D
nuy_global(1:snut3D,1:kchia_q+1)=nuy;%%%%%%%change later
unit_vt_1D=zeros(snut3D,kchia_q+1);
unit_vt_y_1D=zeros(snut3D,kchia_q+1);

PHUONGCHIEUCAO=zeros(snut3D,kchia_q);
for i=1:snut3D
    if i>snut
    PHUONGCHIEUCAO(i,1:index_vertical+1)=H3D(i)/index_vertical*(0:index_vertical);% phuong chieu cao theo do cao 3D
    nuy_local=interp1(H3D_global(i)/kchia_q*(0:kchia_q),nuy_global(i,1:kchia_q+1),PHUONGCHIEUCAO(i,1:index_vertical+1),'linear','extrap'); %noi suy 
    for k=1:index_vertical %chi so chieu cao
        for j=1:2
            [xi,W]=tpso_Gauss_1d(2,j);
            [N1,N2]=shapefunc2D(xi);
            N=[N1 N2];
            L=H3D(i)/index_vertical*h2;%det J
            bien_y_i=N*[PHUONGCHIEUCAO(i,k);PHUONGCHIEUCAO(i,k+1)]*h2;
            
            nuy_1D_i=N*[nuy_local(k);nuy_local(k+1)];
            
            unit_vt_1D_i=N*L/2*W(j);
            unit_vt_1D(i,k)=unit_vt_1D(i,k)+unit_vt_1D_i(1);
            unit_vt_1D(i,k+1)=unit_vt_1D(i,k+1)+unit_vt_1D_i(2);
            
            unit_vt_y_1D_i=N*bien_y_i*L/2*W(j);
            unit_vt_y_1D(i,k)=unit_vt_y_1D(i,k)+unit_vt_y_1D_i(1);
            unit_vt_y_1D(i,k+1)=unit_vt_y_1D(i,k+1)+unit_vt_y_1D_i(2);
        end
    end
    end
end
%% Tinh cac F
F0=zeros(snut3D,1);
F1=zeros(snut3D,1);

for i =1:snut3D
        F0(i)=nuy_global(i,:)*unit_vt_1D(i,:)';
        F1(i)=nuy_global(i,:)*unit_vt_y_1D(i,:)';
end
%% tinh van toc
u1=zeros(snut3D,1);
for i=1:kchia_n+1
    for j=1:kchia_q+1
        u1(luoi3D(1:kchia_m+1,i,j))=omega*(r1+(0:kchia_m)/kchia_m*(r2-r1));
    end
end
A1=zeros(snut3D,1);
B1=zeros(snut3D,1);
B2=zeros(snut3D,1);
for i=1:snut3D
    A1(i)=F1(i)-F1(chiso_max(i))/F0(chiso_max(i))*F0(i);
    B1(i)=u1(i)/F0(chiso_max(i))*F0(i);
    B2(i)=0;
end
vel_u=zeros(snut3D,1);
vel_w=zeros(snut3D,1);
for i=1:snut3D
    vel=[dp_dphi(chiso_min(i)); dp_dr(chiso_min(i))].*A1(i)+[B1(i) B2(i)];
    vel_u(i)=vel(1);
    vel_w(i)=vel(2);
    vel_v(i)=0;%%%%%%%%%%%%%%%%change later
end
% for i=1:snut3D
%     if i>snut
%     PHUONGCHIEUCAO(i,1:index_vertical+1)=H3D(i)/index_vertical*(0:index_vertical);% phuong chieu cao theo do cao 3D
%     nuy_local=interp1(H3D_global(i)/kchia_q*(0:kchia_q),nuy_global(i,1:kchia_q+1),PHUONGCHIEUCAO(i,1:index_vertical+1),'linear','extrap'); %noi suy 
%     for k=1:index_vertical %chi so chieu cao
%         for j=1:2
%             [xi,W]=tpso_Gauss_1d(2,j);
%             [N1,N2]=shapefunc2D(xi);
%             N=[N1 N2];
%             L=H3D(i)/index_vertical*h2;%det J
%             bien_y_i=N*[PHUONGCHIEUCAO(i,k);PHUONGCHIEUCAO(i,k+1)]*h2;
%             
%             nuy_1D_i=N*[nuy_local(k);nuy_local(k+1)];
%             
%             unit_vt_1D_i=N*L/2*W(j);
%             unit_vt_1D(i,k)=unit_vt_1D(i,k)+unit_vt_1D_i(1);
%             unit_vt_1D(i,k+1)=unit_vt_1D(i,k+1)+unit_vt_1D_i(2);
%             
%         end
%     end
%     end
% end
%% tinh nhiet
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
Coordinate3D_polar=zeros(3,snut3D);
Coordinate3D_nature=zeros(3,snut3D);
for i=1:snut3D
    [R,S,T]=ind2sub(size(luoi3D),find(luoi3D==i));
    Coordinate3D_polar(1,i)=pi/2-(S-1)*beta/kchia_n;
    Coordinate3D_polar(2,i)=(r2-(r2-r1)/kchia_m*(R-1));
    Coordinate3D_polar(3,i)=(T-1)/kchia_q*H3D(i);
end

Coordinate3D_nature(1,:)=Coordinate3D_polar(2,:).*cos(Coordinate3D_polar(1,:));
Coordinate3D_nature(2,:)=Coordinate3D_polar(2,:).*sin(Coordinate3D_polar(1,:));
Coordinate3D_nature(3,:)=Coordinate3D_polar(3,:);
node_index3D=zeros(8,spt3D);
for i=1:spt3D
    
    [R,S,T]=ind2sub(size(luoi_phantu3D),find(luoi_phantu3D==i));
    
    node_index3D(1:8,i)=[luoi3D((R+1),S,T),luoi3D((R+1),S+1,T),luoi3D(R,S+1,T),luoi3D(R,S,T),...
            luoi3D((R+1),S,T+1),luoi3D((R+1),S+1,T+1),luoi3D(R,S+1,T+1),luoi3D(R,S,T+1)];
    
end
%%
for i=1:spt3D
    xP=Coordinate3D_nature(1,node_index3D(1:8,i));
    yP=Coordinate3D_nature(3,node_index3D(1:8,i));
    zP=Coordinate3D_nature(2,node_index3D(1:8,i));
    for j=1:8
        [R,S,T,W]=tpso_Gauss_3d(2,j);
        [N1,N2,N3,N4,N5,N6,N7,N8]=hamdang3D(R,S,T);
        N=[N1,N2,N3,N4,N5,N6,N7,N8]';
        dNrst=[-(1/8 - 1/8*S)*(1 - T), -(1/8 - 1/8*R)*(1 - T), -(1/8 - 1/8*R)*(1 - S);...
            [1/8*(1 - S)*(1 - T), -(1/8 + 1/8*R)*(1 - T), -(1/8 + 1/8*R)*(1 - S)];...
            [1/8*(1 + S)*(1 - T), 1/8*(1 + R)*(1 - T), -(1/8 + 1/8*R)*(1 + S)];...
            [-(1/8 + 1/8*S)*(1 - T), 1/8*(1 - R)*(1 - T), -(1/8 - 1/8*R)*(1 + S)];...
            [-(1/8 - 1/8*S)*(1 + T), -(1/8 - 1/8*R)*(1 + T), 1/8*(1 - R)*(1 - S)];...
            [1/8*(1 - S)*(1 + T), -(1/8 + 1/8*R)*(1 + T), 1/8*(1 + R)*(1 - S)];...
            [1/8*(1 + S)*(1 + T), 1/8*(1 + R)*(1 + T), 1/8*(1 + R)*(1 + S)];...
            -(1/8 + 1/8*S)*(1 + T), 1/8*(1 - R)*(1 + T), 1/8*(1 - R)*(1 + S)];...
        J3D=Jctugiac3D(xP,yP,zP,R,S,T);
        Jn=inv(J3D);
        dNx(1:8)=Jn(1,1)*dNrst(1:8,1)+Jn(1,2)*dNrst(1:8,2)+Jn(1,3)*dNrst(1:8,3);
        dNy(1:8)=Jn(2,1)*dNrst(1:8,1)+Jn(2,2)*dNrst(1:8,2)+Jn(2,3)*dNrst(1:8,3);
        dNz(1:8)=Jn(3,1)*dNrst(1:8,1)+Jn(3,2)*dNrst(1:8,2)+Jn(3,3)*dNrst(1:8,3);
        
        A=N'*xP';
        B=N'*yP';
        C=N'*zP';
        U=N'*vel_u(node_index3D(1:8,i));
        V=N'*vel_v(node_index3D(1:8,i))';
        W=N'*vel_w(node_index3D(1:8,i));
%         K=N'*(U*dNx+1/
    end
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

% ve hinh

% 
mtp=zeros(kchia_m+1,kchia_n+1);
% 
mtp(:,:)=vtp(luoi(:,:));
% 
% % [X,Y,Z]=meshgrid(pp(1,:),pp(2,:),H);

% s=surf(X,Y,Z);
figure(2)
plot(pp(1,:),pp(2,:),'+')
figure(1)
x=pp(1,:)';y=pp(2,:)';z=vtp'*h2^2/(12*pi*1500/60*0.039*r2^2);
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
plot(1:kchia_m+1,mtp(kchia_m/2+1,1:kchia_n+1)*h2^2/(12*pi*1500/60*0.039*r2^2))
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


