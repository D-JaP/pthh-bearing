clc
clear all
thu=1;
while (thu<=1)
    if (thu~=1)
clearvars -except u_new thu
    end
m=16;n=16;%so khoang chia
snut=(m+1)*(n+1);
spt=m*n;
phi=linspace(0,2*pi,n+1); % goc phi(phi,k,
Pa=0; % ap suat cap dau (Pa or N)
R=0.035;%(m)
L=0.08;%(m)
nuy=0.008;%do nhot dong luc hoc (N/m^2.s)
omega=300/60*2*pi;%toc do goc (rad/s)
c=0.0003;%khe ho ban kinh (m)
epsilon=0.33;%do lech tam tuong doi
Wy=0;Wx=0.26;%Wx=Wx/(6*nuy*omega*(R/c)^2*R*L);

if (thu==1)
    
x=0.5;y=0.5;u=[x;y];

else
    u=u_new;
    x=u_new(1);
    y=u_new(2);
end
h=1+x*cos(phi)+y*sin(phi);
B=x*sin(phi)-y*cos(phi);
e=epsilon*c;%do lech tam

lamda=0.7;
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
for i=1:(m+1)*(n+1)
    [row,col]=find(luoi==i);
    p(1,i)=(1-col)/n*2*pi;
    p(2,i)=(row-1)/m;
end

t=zeros(4,spt);
for i=1:m*n
        [row,col]=find(luoi_phantu==i);
        t(1:4,i)=[luoi((row+1),col),luoi((row+1),col+1),luoi(row,col+1),luoi(row,col)];
end

%mang chua h cua moi nut
H=zeros(snut,1);
chiso_h=zeros(snut,1);
for i=1:m+1
    for j=1:n+1
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
for i=1:m*n
    xP=p(1,t(1:4,i));
    yP=p(2,t(1:4,i));
    for j=1:4
        switch j
            case 1
                eta=-1/sqrt(3);xi=-1/sqrt(3);
            case 2
                eta=-1/sqrt(3);xi=1/sqrt(3);
            case 3
                eta=1/sqrt(3);xi=1/sqrt(3);
            case 4
                eta=1/sqrt(3);xi=-1/sqrt(3);
        end
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
    
    HH3=(N'*H(t(1:4,i))).^3;
    
    
    h2cos=N'*(H(t(1:4,i)).^2.*cos(phi(chiso_h(t(1:4,i))))');
    h2sin=N'*(H(t(1:4,i)).^2.*sin(phi(chiso_h(t(1:4,i))))');
    %tinh cac ma tran do cung 
    Ai=HH3*(dNx'*dNx+dNz'*dNz)*abs(det(J));
    Adx=h2cos*(dNx'*dNx+dNz'*dNz)*abs(det(J));
    Ady=h2sin*(dNx'*dNx+dNz'*dNz)*abs(det(J));

    mtA(t(1:4,i),t(1:4,i))=mtA(t(1:4,i),t(1:4,i))+Ai;
    mtAdx(t(1:4,i),t(1:4,i))=mtAdx(t(1:4,i),t(1:4,i))+Adx;
    mtAdy(t(1:4,i),t(1:4,i))=mtAdy(t(1:4,i),t(1:4,i))+Ady;
    
    %tinh vecto tai
    Bj=N*N'*B(chiso_h(t(1:4,i)))'*abs(det(J));
    Bx=N*(N'*sin(phi(chiso_h(t(1:4,i))))')*abs(det(J));
    By=N*(-N'*cos(phi(chiso_h(t(1:4,i))))')*abs(det(J));
    
    b(t(1:4,i))=b(t(1:4,i))+Bj;
    bnx(t(1:4,i))=bnx(t(1:4,i))+Bx;
    bny(t(1:4,i))=bny(t(1:4,i))+By;
    
    %tinh vecto S R
    SS=N*(N'*cos(phi(chiso_h(t(1:4,i))))')*abs(det(J));
    RR=N*(N'*sin(phi(chiso_h(t(1:4,i))))')*abs(det(J));

    vt_S(t(1:4,i))=vt_S(t(1:4,i)) + SS;
    vt_R(t(1:4,i))=vt_R(t(1:4,i)) + RR;
    end
end

%tao vec to ap suat 
vtp=zeros(snut,1);
q=1;


for i=1:m+1
    for j=1:n+1
        if (i~=1&&i~=m+1)
            vectoan(q)=luoi(i,j);
            q=q+1;
            vectoan=sort(vectoan);
        end
    end
end

Ast=mtA(vectoan,vectoan);
bst=b(vectoan);
vtp(vectoan)=inv(Ast)*bst;
cd=linspace(0,1,m+1);
mtp(:,:)=vtp(luoi(:,:));


demIa=1;
demIb=1;

for i=1:m+1
    for j=1:n+1
        if (vtp(luoi(i,j))>0&&i~=1&&i~=m+1)
            Ia(demIa)=luoi(i,j);
            demIa=demIa+1;
            Ia=sort(Ia);
        end
        
        if (vtp(luoi(i,j))<=0)
            Ib(demIb)=luoi(i,j);
            vtp(luoi(i,j))=0;
            mtp(i,j)=0;
            demIb=demIb+1;
            Ib=sort(Ib);
        end
    end
end

Ast=mtA(Ia,Ia);
bst=b(Ia);
vtp(:)=0;
vtp(Ia)=Ast\bst;
mtp(:,:)=0;
mtp(:,:)=vtp(luoi(:,:));
%ve do thi
[X,Y]=meshgrid(cd,phi);
figure(1)
surf(X',Y',mtp*6*nuy*omega*(R/c)^2/1000)

% tim cac nut thuoc ia ib
% 
% for k=0:(m+1)*(n+1);
%     
%     if (k==0)
%         Ia=1:(m+1)*(n+1);
%         Ib=0;
%     end
% 
% Am=zeros((m+1)*(n+1),(m+1)*(n+1));
% %xac dinh A mu
% for i=1:(m+1)*(n+1)
%     for j=1:(m+1)*(n+1)
%         if (sum(Ia==i)==1&&sum(Ia==j)==1)
%         Am(i,j)=mtA(i,j);
%         end
%         if (i==j)
%             Am(i,j)=1;
%         end
%     end
% end
% % xac dinh b mu
% Bm=zeros((m+1)*(n+1),1);
% for i=1:(m+1)*(n+1)
%     if (sum(Ia==i)==1)
%         Bm(i)=b(i);
%     end
%     if (sum(Ib==i)==1)
%         Bm(i)=0;
%     end
% end
% 
% p_nga=inv(Am)*Bm;
% q_nga=mtA*p_nga-b;
% if (sum((p_nga(:)>=0)==0)==1&&sum((q_nga(:)>=0)==0))
%     k
%     break
%     
% else
%     demIa=1;demIb=1;
%     for i=1:(m+1)*(n+1)
%         if (p_nga(i)>0||q_nga(i)<0)
%             Ia(demIa)=i;
%             demIa=demIa+1;
%         else
%             Ib(demIb)=i;
%             demIb=demIb+1;
%         end
%     end
% end
% end



% chuyen ma tran A-> Amu, b->bm
Am=zeros(snut,snut);
bm=zeros(snut,1);
for i=1:snut
    for j=1:snut
        if (sum(sum(i==Ia))==1&&sum(sum(j==Ia))==1)
            Am(i,j)=mtA(i,j);
        end
        if (sum(sum(i==Ib))==1)
            Am(i,i)=1;
        end    
        if ((sum(sum(i==Ia))==1&&sum(sum(j==Ib))==1)||(sum(sum(i==Ib))==1&&sum(sum(j==Ia))==1))
            Am(i,j)=0;         
        end     
    end
end
%tim bmu
for i=1:snut
    if (sum(sum(i==Ia))==1)
        bm(i)=b(i);
    end
    if (sum(sum(i==Ib))==1)
        bm(i)=0;
    end
end

vtp=Am\bm;

mtAmx=zeros(snut,snut);
mtAmy=zeros(snut,snut);
bmx=zeros(snut,1);
bmy=zeros(snut,1);
for i=1:snut
    for j=1:snut
        if (sum(sum(i==Ia))==1&&sum(sum(j==Ia))==1)
            mtAmx(i,j)=mtAdx(i,j);
            mtAmy(i,j)=mtAdy(i,j);
        else 
            mtAmx(i,j)=0;
            mtAmy(i,j)=0;
            bmx(i)=0;
            bmy(i)=0;
        end
    end
end

for i=1:snut
    if (sum(sum(i==Ia))==1)
        bmx(i)=bnx(i);
        bmy(i)=bny(i);
    end
end
pxy=Am\[-mtAmx*vtp+bmx -mtAmy*vtp+bmy];
% for i=1:(m)*(n)
%     for j=1:4
%         
% mtpnx(i,j)=pxy(t(j,i),1);
% mtpny(i,j)=pxy(t(j,i),2);
%     end
% end
%D_jc=-[S'*pnx' S'*pny';R'*pnx' R'*pnx'];
D_jc=-[vt_S'*pxy(:,1) vt_S'*pxy(:,2);...
       vt_R'*pxy(:,1) vt_R'*pxy(:,2)];

mtfS=-vt_S'*vtp;
mtfR=-vt_R'*vtp;
  
fx=sum(mtfS);
fy=sum(mtfR);


u_new=u-D_jc\([fx;fy]-[Wx;Wy])
thu=thu+1
if ((norm(u_new-u)/norm(u_new))<=1e-5&&(norm([fx;fy]-[Wx;Wy])/norm([Wx;Wy]))<1e-5)
    disp('done')
break
end
end
