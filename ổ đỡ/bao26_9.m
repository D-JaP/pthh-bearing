clc
clear all
thu=1;
Pc=123*10^3;
m=30;n=80;%so khoang chia
snut=(m+1)*(n);
spt=m*n;
phi=linspace(0,2*pi,n+1); % goc phi(phi,k,
Pa=0; % ap suat cap dau (Pa or N)
R=0.05;%(m)
L=0.05;%(m)
nuy=0.015;%do n0hot dong luc hoc (N/m^2.s)
omega=300/60*(2*pi);%toc do goc (rad/s)
c=0.00005;%khe ho ban kinh (m)
Rcd=0.005;
goc_cd=asin(Rcd/R);
lamda=L/(2*R);
% danh so nut
p_so=1;
luoi=zeros(m+1,n);
Pcd=10000/(6*nuy*omega*(R/c)^2);
for j=1:n
    if mod(j,2)~=0
        for i=1:m+1
            luoi(i,j)=p_so;
            p_so=p_so+1;
        end
    end
    if mod(j,2)==0 
        for i=1:m+1
            luoi(i,j)=p_so;
            p_so=p_so+1;
        end
    end
end
%danh so phan tu
p_so2=1;
luoi_phantu=zeros(m,n);
for j=1:n
%     if mod(i,2)~=0
        for i=1:m
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
% xac dinh cac nut thuoc Ia Ib
%         Ia=luoi(1:m+1,1:n/2+1);
%         Ib=luoi(1:m+1,n/2+2:n+1);
%==========================================================================================
p=zeros(2,snut);
for i=1:snut
    [row,col]=find(luoi==i);
    p(1,i)=(1-col)/n*2*pi;
    p(2,i)=(1/2-(m-row+1)/m)*L/(2*R)*2;
end
t=zeros(4,spt);
for i=1:m*n
        [row,col]=find(luoi_phantu==i);
        if (col~=n)
            t(1:4,i)=[luoi((row+1),col),luoi((row+1),col+1),luoi(row,col+1),luoi(row,col)];
        else 
            t(1:4,i)=[luoi((row+1),col),luoi((row+1),1),luoi(row,1),luoi(row,col)];
        end
end







while (1)


Wy=0;
% Wx=-(140)/(6*nuy*omega*(R/c)^2*R*R);
% Stormfield=(nuy*omega/2/pi*L*2*R*(R/c)^2)/(70+70);
% Stormfield=0.2;
% Wx=(L/2/R*1/3/pi/Stormfield)
Stormfield=2;
Wx=(L/2/R*1/3/pi/Stormfield)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (thu==1)

    x=0.5;
    y=0.5;
    u=[x;y];

    else

        u=u_new;
        x=u_new(1);
        y=u_new(2);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1+x*cos(phi)+y*sin(phi);
klr=860;%khoi luong rieng
Cp=2000;

%mang chua h cua moi nut
H=zeros(snut,1);
chiso_h=zeros(snut,1);

% tao mang chieu day mang dau
    for i=1:m+1
        for j=1:n
            H(luoi(i,j))=h(j);
            chiso_h(luoi(i,j))=j;
        end
    end



% tich phan so 
mtA=zeros(snut,snut);
mtAdx=zeros(snut,snut);
mtAdy=zeros(snut,snut);


b=zeros(snut,1);
bny=zeros(snut,1);
bnx=zeros(snut,1);
vt_S=zeros(snut,1);
vt_R=zeros(snut,1);

% dieu kien bien khoi tao
q=1;
dem_Ib=1;
dem_Ic=1;
dem_Icd=1;
% for j=1:n
%     if (phi(j)>235/180*pi)
%         bien_capdau=j-1;
%     break
%     end
% end
i=1;
    for j=1:n
        if ((phi(j)>=80/180*pi&&phi(j)<=100/180*pi)||(phi(j)>=260/180*pi&&phi(j)<=280/180*pi))
            bien_capdau(i)=j;
            i=i+1;
        end
    end
    bien_cd=(luoi(:,bien_capdau));
    
q=1;
    for i=1:m+1
        for j=1:n
            if (i==1||i==m+1||sum(j==bien_capdau)==1)
                I_boundry(dem_Ib)=luoi(i,j);
                dem_Ib=dem_Ib+1;
                
            else
                I_active(q)=luoi(i,j);
                q=q+1;
 
    %         else 
    %             I_cavitation(dem_Ic)=luoi(i,j);
    %             dem_Ic=dem_Ic+1;
            end

        end
    end
I_cavitation=double.empty(0,0);
%tich phan so Gauss

    xP=p(1,t(1:4,1));
    yP=p(2,t(1:4,1));
    for i=1:spt
    %     xP=p(1,t(1:4,i));
    %     yP=p(2,t(1:4,i));
        if i==spt/2
           chiso_h(luoi(:,1))=n+1;
        end
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

        HH3=(N'*H(t(1:4,i)))^3;
        h2cos=(N'*H(t(1:4,i)))^2*(cos(N'*phi(chiso_h(t(1:4,i)))'));
        h2sin=(N'*H(t(1:4,i)))^2*(sin(N'*phi(chiso_h(t(1:4,i)))'));
        %tinh cac ma tran do cung
        Ai=HH3*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
        Adx=3*h2cos*(dNx'*dNx+dNz'*dNz)*abs(det(J))*W(j);
        Ady=3*h2sin*(dNx'*dNx+dNz'*dNz)*abs(det(J))*W(j);

        mtA(t(1:4,i),t(1:4,i))=mtA(t(1:4,i),t(1:4,i))+Ai;
        mtAdx(t(1:4,i),t(1:4,i))=mtAdx(t(1:4,i),t(1:4,i))+Adx;
        mtAdy(t(1:4,i),t(1:4,i))=mtAdy(t(1:4,i),t(1:4,i))+Ady;

        %tinh vecto tai
        %Bj=12*pi*N*(N'*dH(t(1:4,i))')*abs(det(J))*W(j);
        Bj=N*(x*sin(N'*phi(chiso_h(t(1:4,i)))')-y*cos(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J))*W(j);
        %Bj=12*pi*N*(dNx*H(t(1:4,i)))*abs(det(J))*W(j);
        Bx=N*(sin(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J))*W(j);
        By=-N*(cos(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J))*W(j);

        b(t(1:4,i))=b(t(1:4,i))+Bj;
        bnx(t(1:4,i))=bnx(t(1:4,i))+Bx;
        bny(t(1:4,i))=bny(t(1:4,i))+By;
    %     
    %     
        SS=N*(cos(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J))*W(j);
        RR=N*(sin(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J))*W(j);

        vt_S(t(1:4,i))=vt_S(t(1:4,i)) + SS;
        vt_R(t(1:4,i))=vt_R(t(1:4,i)) + RR;


        end
    end
chiso_h(luoi(:,1))=1;

 

    while (1)

    % % thiet lap Am
    Am=zeros(snut,snut);
    Bm=zeros(snut,1);


    Am(I_active,I_active)=mtA(I_active,I_active);
    Bm(I_active)=b(I_active);

    
    for i=1:snut
        if (sum(I_active==i)==0) % i khong thuoc I active
            Am(i,i)=1;
        end
    end

    %tinh ap suat

    vtp=Am\Bm;
    
    vtq=zeros(snut,1);

    vtq([I_cavitation,I_active])=mtA([I_cavitation,I_active],([I_cavitation,I_active]))*vtp([I_cavitation,I_active])-b([I_cavitation,I_active]);

    if (sum(vtp(:)<-1e-7)==0&&(sum(vtq(:)<=-1e-7))==0)
            break
    end

    clear I_active
    clear I_cavitation

    dem_Ia=1;

    for i=1:snut
        if ((vtp(i)>0&&abs(vtp(i))>1e-7)||(vtq(i)<0&&abs(vtq(i))>1e-7))
            I_active(dem_Ia)=i;
            dem_Ia=dem_Ia+1;
        end
    end
    dem_Ic=1;
    for i=1:snut
        if (sum(i==I_active)==0&&sum(i==I_boundry)==0)
            I_cavitation(dem_Ic)=i;
            dem_Ic=dem_Ic+1;
        end
    end

    end
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

Am_X=zeros(snut,snut);
Am_Y=zeros(snut,snut);
Bm_X=zeros(snut,1);
Bm_Y=zeros(snut,1);

Am_X(I_active,I_active)=mtAdx(I_active,I_active);
Am_Y(I_active,I_active)=mtAdy(I_active,I_active);
Bm_X(I_active)=bnx(I_active);
Bm_Y(I_active)=bny(I_active);

% for i=1:snut
%     if (sum(I_active==i)==0) % i khong thuoc I active
%         Am_X(i,i)=1;
%         Am_Y(i,i)=1;
%     end
% end

%tinh dp/dx, dp/dy
pxy=Am\[-Am_X*vtp+Bm_X -Am_Y*vtp+Bm_Y];

pnx=pxy(:,1);
pny=pxy(:,2);

D_jc=-[vt_S'*pnx vt_S'*pny;...
       vt_R'*pnx vt_R'*pny]

fx=-vt_S'*vtp
fy=-vt_R'*vtp

u_new=u-D_jc\([fx;fy]-[Wx;Wy])

if ((norm(u_new-u)/norm(u_new))<=1e-5&&(norm([fx;fy]-[Wx;Wy])/norm([Wx;Wy]))<1e-5)
    disp('done')
break
end

u=u_new;
thu=thu+1;
end
% disp('Pmax:'),max(vtp)
% cd=linspace(-1,1,m+1);
% mtp=zeros(m+1,n+1);
% mtq=zeros(m+1,n+1);
% mtp(:,1:n)=vtp(luoi(:,1:n));
% mtq(:,1:n)=vtq(luoi(:,1:n));
% mtp(:,n+1)=mtp(:,1);
% mtq(:,n+1)=mtq(:,1);
% [X,Y]=meshgrid(cd,phi);
% figure(1)
% surf(X',Y',mtp)
% axis tight
% xlabel('L(m)')
% ylabel('\theta (rad)')
% % zlabel('P(bar)')
% hh = rotate3d;
% set(hh,'ActionPreCallback',...
%     'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(hh,'ActionPostCallback',...
%     'set(gcf,''windowbuttonmotionfcn'','''')')
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
% figure(6)
% plot(phi(1:n),mtp(m/2,[n/2:n,1:n/2-1])*6*nuy*omega*(R/c)^2)
% figure(7) 
% plot(cd,mtp(:,bien_capdau))
% hold on
% plot(cd,mtp(:,bien_capdau+3))
% plot(cd,mtp(:,bien_capdau+2))
% plot(cd,mtp(:,bien_capdau+1))
K=6*nuy*omega*R^2*(R/c)^2*D_jc/(Wx*6*nuy*omega*(R/c)^2*R*R)