clc
clear all
close all


thu=1;
input_o_chan_nhiet
%%
global snut snut3D spt2D spt3D 
snut=(kchia_m+1)*(kchia_n+1);
snut3D=(kchia_m+1)*(kchia_n+1)*(kchia_q+1);
spt2D=kchia_m*kchia_n;
spt3D=kchia_m*kchia_n*kchia_q;
%% danh so nut luoi 2D
global luoi
[luoi]=danhsonutluoi2D(kchia_m,kchia_n);
%% danh so nut luoi 3D
global luoi3D
[luoi3D]=danhsonutluoi3D(kchia_m,kchia_n,kchia_q);% danh so nut

%% danh so phan tu 2D 3D
global luoi_phantu luoi_phantu3D
luoi_phantu=luoiphantu2D(kchia_m,kchia_n);

luoi_phantu3D=luoiphantu3D(kchia_m,kchia_n,kchia_q);

%% toa do va luoi phan tu 2D
global p
p=tinhtoadophantu(snut,luoi,kchia_m,kchia_n,beta,r1,r2);
%% tinh toa do tu nhien
global pp
pp=zeros(2,snut);
pp(1,:)=p(2,:).*cos(p(1,:));
pp(2,:)=p(2,:).*sin(p(1,:));
%% chuyen cac nut vao phan tu
global t
t=chuyennutvaophantu(spt2D,luoi_phantu,kchia_m,kchia_n,luoi);


%%

    %%
    global h1 h2 deltah h2_new delta
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
        
%% tinh h
congthuc_ochan_nhiet
    h=zeros(kchia_m+1,kchia_n+1);
%     phuongchieucao=zeros(kchia_n+1,kchia_n+1,kchia_n+1);%%%%%%%%%%%%%%%%%%
%     PHUONGCHIEUCAO=zeros((kchia_n+1)*(kchia_n+1),kchia_n+1);%%%%%%%%%%%%%%%%%%
    for i=1:kchia_m+1
        for j=1:kchia_n+1
            h(i,j)=eval(h_ct(pp(1,luoi(i,j)),pp(2,luoi(i,j))));
%             h(i,j)=1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/beta);
%             h(i,j)=(1 + delta*(1 - acot(pp(2,luoi(i,j))/pp(1,luoi(i,j)))/sqrt(pp(1,luoi(i,j))^2 + pp(2,luoi(i,j))^2)));
            
%             phuongchieucao(i,j,1:kchia_q+1)=(0:kchia_q)*h(i,j);%%%%%%%%%%%%%%%%
%             PHUONGCHIEUCAO(luoi(i,j),1:kchia_q+1)=(0:kchia_q)*h(i,j);
            
        end
    end
%%
    [I_active,I_boundry]=allocate_element(kchia_m,kchia_n,luoi);    %xac dinh I_active,I_boundary
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
    mtA=zeros(snut,snut);%mtA matran do cung tong quat
    mtA_d=zeros(snut,snut);%matran do cung khi them dieu kien bien
    mtAh=zeros(snut,snut);% ma tran do cung dao ham theo h
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

%% tich phan so Gauss
    
[mtA,mtAh,dA_dphi,dA_dr,db_dphi,db_dr,b,bh,unit_vt,vt_X,vt_Y]=tichphanGauss2D_ochan_nhiet();

    %%
for i=1:snut
    if (sum(I_active==i)==0) % i khong thuoc I active
        mtA_d(i,i)=1;
    end
end
%% them dk bien
mtA_d(I_active,I_active)=mtA(I_active,I_active);
mtAh_d(I_active,I_active)=mtAh(I_active,I_active);
dmtA_dphi(I_active,I_active)=dA_dphi(I_active,I_active);
dmtA_dr(I_active,I_active)=dA_dr(I_active,I_active);
%%
b(I_boundry)=0;
bh(I_boundry)=0;
db_dphi(I_boundry)=0;
db_dr(I_boundry)=0;

vtp=(mtA_d\b)';
Fz=unit_vt'*vtp';

M_X=vt_X'*vtp';
M_Y=vt_Y'*vtp';

Xc_tinh=M_X/Fz
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
%% chuyen doi medium to max
chiso_max=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
chiso_min=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
for k=1:kchia_q+1
    for j=1:kchia_n+1
        for i=1:kchia_m+1
            chiso_max(luoi3D(i,j,k))=luoi3D(i,j,kchia_q+1);
            chiso_min(luoi3D(i,j,k))=luoi3D(i,j,1);
            
        end 
    end
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

%% toa do va luoi phan tu



Coordinate3D_polar=toadophantu3D(snut3D,kchia_m,kchia_n,kchia_q,r1,r2,beta,H3D,luoi3D);
Coordinate3D_nature=zeros(3,snut3D);
Coordinate3D_nature(1,:)=Coordinate3D_polar(2,:).*cos(Coordinate3D_polar(1,:));
Coordinate3D_nature(2,:)=Coordinate3D_polar(2,:).*sin(Coordinate3D_polar(1,:));
Coordinate3D_nature(3,:)=Coordinate3D_polar(3,:);
%%
node_index3D=chuyennutvaophantu3D(spt3D,luoi_phantu3D,luoi3D);
%%
for i=1:spt3D
    xP=Coordinate3D_nature(1,node_index3D(1:8,i));
    yP=Coordinate3D_nature(3,node_index3D(1:8,i));
    zP=Coordinate3D_nature(2,node_index3D(1:8,i));
    for j=1:8
        [R,S,T,W]=tpso_Gauss_3d(2,j);
        [N1,N2,N3,N4,N5,N6,N7,N8]=hamdang3D(R,S,T);
        N=[N1,N2,N3,N4,N5,N6,N7,N8]';
        dNrst=dNrst_3D(R,S,T);
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
%%
vehinh


