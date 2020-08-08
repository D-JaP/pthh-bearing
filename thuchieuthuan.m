clc
clear all
m=100;n=100;%so khoang chia
snut=(m+1)*(n+1);
spt=m*n;
phi=linspace(0,2*pi,n+1); % goc phi(phi,k,
Pa=0; % ap suat cap dau (Pa or N)
R=0.0508;%(m)
L=0.1016;%(m)
nuy=1.3*10^-6;%do n0hot dong luc hoc (N/m^2.s)
omega=300/60*2*pi;%toc do goc (rad/s)
c=5.08*10-5;%khe ho ban kinh (m)
epsilon=0.5;%do lech tam tuong doi
lamda=L/(2*R);


h=1+epsilon*cos(phi);
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
    p(2,i)=(m+1-row)/m;
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
% mtAdx=zeros(snut,snut);
% mtAdy=zeros(snut,snut);
% mtAdX=zeros(snut,snut);
% mtAdY=zeros(snut,snut);

b=zeros(snut,1);
% bny=zeros(snut,1);
% bnx=zeros(snut,1);
% vt_S=zeros(snut,1);
% vt_R=zeros(snut,1);

% dieu kien bien
q=1;
for i=1:m+1 
    for j=1:n+1
        if (i~=1&&i~=m+1&&phi(j)<190/360*2*pi&&j~=1)
            vectoan(q)=luoi(i,j);
            I_active(q)=luoi(i,j);
            q=q+1;
            vectoan=sort(vectoan);
        end
    end
end

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
%     h2cos=(N'*H(t(1:4,i)))^2*(cos(N'*phi(chiso_h(t(1:4,i)))'));
%     h2sin=(N'*H(t(1:4,i)))^2*(sin(N'*phi(chiso_h(t(1:4,i)))'));
    %tinh cac ma tran do cung
    Ai=HH3*(dNx'*dNx+1/(4*lamda^2)*(dNz'*dNz))*abs(det(J));
%     Adx=h2cos*(dNx'*dNx+dNz'*dNz)*abs(det(J));
%     Ady=h2sin*(dNx'*dNx+dNz'*dNz)*abs(det(J));

    mtA(t(1:4,i),t(1:4,i))=mtA(t(1:4,i),t(1:4,i))+Ai;
%     mtAdx(t(1:4,i),t(1:4,i))=mtAdx(t(1:4,i),t(1:4,i))+Adx;
%     mtAdy(t(1:4,i),t(1:4,i))=mtAdy(t(1:4,i),t(1:4,i))+Ady;
    
    %tinh vecto tai
    %Bj=12*pi*N*(N'*dH(t(1:4,i))')*abs(det(J));
    Bj=12*pi*N*(epsilon*sin(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J));
    %Bj=12*pi*N*(dNx*H(t(1:4,i)))*abs(det(J));
%     Bx=N*(sin(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J));
%     By=-N*(cos(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J));
    
    b(t(1:4,i))=b(t(1:4,i))+Bj;
%     bnx(t(1:4,i))=bnx(t(1:4,i))+Bx;
%     bny(t(1:4,i))=bny(t(1:4,i))+By;
%     
%     
%     SS=N*(cos(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J));
%     RR=N*(sin(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J));

%     vt_S(t(1:4,i))=vt_S(t(1:4,i)) + SS;
%     vt_R(t(1:4,i))=vt_R(t(1:4,i)) + RR;
    end
end


%tao vec to ap suat 
vtp=zeros(snut,1);





% tinh ap suat 
Ast=mtA(vectoan,vectoan);
bst=b(vectoan);
vtp(vectoan)=Ast\bst;
cd=linspace(0,1,m+1);
mtp(:,:)=vtp(luoi(:,:));

[X,Y]=meshgrid(cd,phi);
figure(1)
surf(X',Y',mtp)

disp('Pmax:'),max(vtp)
vtp_giua=vtp(luoi(m/2,:));
figure(2)
plot(phi,vtp_giua)
phi(find(vtp_giua==max(vtp_giua)));