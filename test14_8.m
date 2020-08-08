clc
clear all
thu=1;
while (1)
    if (thu~=1)
clearvars -except u_new thu
    end
m=50;n=100;%so khoang chia
snut=(m+1)*(n+1);
spt=m*n;
phi=linspace(0,2*pi,n+1); % goc phi(phi,k,
Pa=0; % ap suat cap dau (Pa or N)
R=0.04;%(m)
L=0.08;%(m)
nuy=0.008;%do n0hot dong luc hoc (N/m^2.s)
omega=300/60*2*pi;%toc do goc (rad/s)
c=0.0003;%khe ho ban kinh (m)
epsilon=0.5;%do lech tam tuong doi
lamda=L/(2*R);


Wy=0;
Wx=(70+70)/(6*nuy*omega*(R/c)^2);
%Wx=0.573;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (thu==1)
    
x=0.3;
y=0;
u=[x;y];

else
    u=u_new;
    x=u_new(1);
    y=u_new(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1+x*cos(phi)+y*sin(phi);
dh=epsilon*sin(phi);
e=epsilon*c;%do lech tam


klr=860;%khoi luong rieng
Cp=2000;
% danh so nut
p_so=1;
luoi=zeros(m+1,n+1);
for j=1:n+1
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
        t(1:4,i)=[luoi((row+1),col),luoi((row+1),col+1),luoi(row,col+1),luoi(row,col)];
end

%mang chua h cua moi nut
H=zeros(snut,1);
chiso_h=zeros(snut,1);

% tao mang chieu day mang dau
for i=1:m+1
    for j=1:n+1
        H(luoi(i,j))=h(j);
        dH(luoi(i,j))=dh(j);
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
for i=1:m+1 
    for j=1:n+1
        if (i~=1&&i~=m+1&&j~=1&&j~=n+1)
            
            I_active(q)=luoi(i,j);
            q=q+1;
            
        elseif (i==1||i==m+1||j~=1||j~=n+1)
            I_boundry(dem_Ib)=luoi(i,j);
            dem_Ib=dem_Ib+1;
%         else 
%             I_cavitation(dem_Ic)=luoi(i,j);
%             dem_Ic=dem_Ic+1;
        end
        
    end
end
I_cavitation=double.empty(0,0);
%tich phan so Gauss
for i=1:spt
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
    
    HH3=(N'*H(t(1:4,i)))^3;
    h2cos=(N'*H(t(1:4,i)))^2*(cos(N'*phi(chiso_h(t(1:4,i)))'));
    h2sin=(N'*H(t(1:4,i)))^2*(sin(N'*phi(chiso_h(t(1:4,i)))'));
    %tinh cac ma tran do cung
    Ai=HH3*(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
    Adx=h2cos*(dNx'*dNx+dNz'*dNz)*abs(det(J))*W(j);
    Ady=h2sin*(dNx'*dNx+dNz'*dNz)*abs(det(J))*W(j);

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


 

while (1)
    
% % thiet lap Am
Am=zeros(snut,snut);
Bm=zeros(snut,1);


Am(I_active,I_active)=mtA(I_active,I_active);
Bm(I_active)=b(I_active);
Am(I_active,[I_cavitation,I_boundry])=0;
Am([I_cavitation,I_boundry],I_active)=0;

for i=1:snut
    if (sum(I_active==i)==0) % i khong thuoc I active
        Am(i,i)=1;
    end
end

%tinh ap suat
disp('aaaaaaaaaaaa')
vtp=Am\Bm;
vtq=zeros(snut,1);

vtq([I_cavitation,I_active])=mtA([I_cavitation,I_active],([I_cavitation,I_active]))*vtp([I_cavitation,I_active])-b([I_cavitation,I_active]);
disp('bbbbbbbbbbbb')
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


disp('Pmax:'),max(vtp)
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
Am_X(I_active,[I_cavitation,I_boundry])=0;
Am_X([I_cavitation,I_boundry],I_active)=0;
Am_Y(I_active,[I_cavitation,I_boundry])=0;
Am_Y([I_cavitation,I_boundry],I_active)=0;
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
       vt_R'*pnx vt_R'*pny];

fx=-vt_S'*vtp
fy=-vt_R'*vtp;

u_new=u-D_jc\([fx;fy]-[Wx;Wy])

if ((norm(u_new-u)/norm(u_new))<=1e-5&&(norm([fx;fy]-[Wx;Wy])/norm([Wx;Wy]))<1e-5)
    disp('done')
break
end
thu=thu+1;
end

cd=linspace(0,1,m+1);
mtp(:,:)=vtp(luoi(:,:));
mtq(:,:)=vtq(luoi(:,:));
[X,Y]=meshgrid(cd,phi);
figure(1)
surf(X',Y',mtp*6*nuy*omega*(R/c)^2)
figure(2)
surf(X',Y',mtq*6*nuy*omega*(R/c)^2)
[row,col]=find(mtp==max(vtp));
phi(col)/pi*180
