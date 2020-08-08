clc
clear all

thu=1;
Pc=123*10^3;
m=40;n=160;%so khoang chia
snut=(m+1)*(n);
spt=m*n;
phi=linspace(0,2*pi,n+1); % goc phi(phi,k,
Pa=0; % ap suat cap dau (Pa or N)
R=0.035;%(m)
L=0.05;%(m)
nuy=0.015;%do n0hot dong luc hoc (N/m^2.s)
omega=300/60*2*pi;%toc do goc (rad/s)
c=0.00005;%khe ho ban kinh (m)
Rcd=0.0025;
goc_cd=asin(Rcd/R);
lamda=L/(2*R);
% danh so nut
p_so=1;
luoi=zeros(m+1,n);
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

% for lap_tong=[0.1:0.1:1,2:1:10]
%     if lap_tong<=1
%     ll=int8(lap_tong*10);
%     else 
%         ll=10+lap_tong-1;
%     end
    Wy=(0)/(6*nuy*omega*(R/c)^2*R*R)
    Wx=(-140)/(6*nuy*omega*(R/c)^2*R*R)
%     Stormfield=(nuy*omega/2/pi*L*2*R*(R/c)^2)/(320)
%     Stormfield=lap_tong;
%     Wx=(L/2/R*1/3/pi/Stormfield);
%     Wx=0.26525823;




while (1)
    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (thu==1)
       
        x=-0.06;
        y=-0.25;
        u=[x;y];
        
    else
        u=u_new;
        x=u_new(1);
        y=u_new(2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    bnx_dot=zeros(snut,1);
    bny_dot=zeros(snut,1);
    vt_S=zeros(snut,1);
    vt_R=zeros(snut,1);
    
    % dieu kien bien khoi tao
    q=1;
    dem_Ib=1;
    dem_Ic=1;
    % for j=1:n
    %     if (phi(j)>235/180*pi)
    %         bien_capdau=j-1;
    %     break
    %     end
    % end
%     for j=1:n
%         if (phi(j)>=225/180*pi)
%             bien_capdau=j;
%             break
%         end
%     end
%     
%     
%% lo cap dau

    q=1;bien_capdau=(180+45)/360*n;
    bien_cd=[];
    for i=1:m+1
        for j=1:n
            if (phi(j)>=phi(bien_capdau)-goc_cd&&phi(j)<=phi(bien_capdau)+goc_cd&&((p(2,luoi(i,j))-p(2,luoi(m/2+1,1)))*R)^2+((p(1,luoi(i,j))-p(1,luoi(m/2,bien_capdau)))*R)^2<=Rcd^2)
                bien_cd(q)=luoi(i,j);

                q=q+1;
            end
        end
    end
    %%
% bien_cd=[];
%     for j=1:n
%         
%         if (((phi(j)/pi*180)>80&&(phi(j)/pi*180)<100)||((phi(j)/pi*180)>260&&(phi(j)/pi*180)<280))
%             bien_cd(q)=j;
%             q=q+1;
%         end
%     end
%     
    q=1;
    I_active=[];
    I_boundry=[];
    for i=1:m+1
        for j=1:n
            
            if (i~=1&&i~=m+1&&sum(luoi(i,j)==bien_cd)==0)
                
                I_active(q)=luoi(i,j);
                q=q+1;
                
            elseif (i==1||i==m+1||sum(luoi(i,j)==bien_cd)~=0)
                I_boundry(dem_Ib)=luoi(i,j);
                
                if (sum(luoi(i,j)==bien_cd)~=0)
                   disp('chan') 
                    
                end
                dem_Ib=dem_Ib+1;
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
            Bx_dot=-2*N*(cos(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J))*W(j);
            By_dot=-2*N*(sin(N'*phi(chiso_h(t(1:4,i)))'))*abs(det(J))*W(j);
            
            b(t(1:4,i))=b(t(1:4,i))+Bj;
            bnx(t(1:4,i))=bnx(t(1:4,i))+Bx;
            bny(t(1:4,i))=bny(t(1:4,i))+By;
            bnx_dot(t(1:4,i))=bnx(t(1:4,i))+Bx_dot;
            bny_dot(t(1:4,i))=bny(t(1:4,i))+By_dot;
            
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
        Am(I_active,[I_cavitation,I_boundry])=0;
        Am([I_cavitation,I_boundry],I_active)=0;
        
        for i=1:snut
            if (sum(I_active==i)==0) % i khong thuoc I active
                Am(i,i)=1;
                if (sum(bien_cd==i)~=0)
                   Bm(i)=0.05 ;
                end
            end
        end
        
        
        %tinh ap suat
        
        vtp=Am\Bm;
        
        vtq=zeros(snut,1);
        
        vtq([I_cavitation,I_active])=(mtA([I_cavitation,I_active],([I_cavitation,I_active]))*vtp([I_cavitation,I_active])-b([I_cavitation,I_active]));
        
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
    Bm_X_dot=zeros(snut,1);
    Bm_Y_dot=zeros(snut,1);
    
    Am_X(I_active,I_active)=mtAdx(I_active,I_active);
    Am_Y(I_active,I_active)=mtAdy(I_active,I_active);
    Am_X(I_active,[I_cavitation,I_boundry])=0;
    Am_X([I_cavitation,I_boundry],I_active)=0;
    Am_Y(I_active,[I_cavitation,I_boundry])=0;
    Am_Y([I_cavitation,I_boundry],I_active)=0;
    Bm_X(I_active)=bnx(I_active);
    Bm_Y(I_active)=bny(I_active);
    Bm_X_dot(I_active)=bnx_dot(I_active);
    Bm_Y_dot(I_active)=bny_dot(I_active);
    
    % for i=1:snut
    %     if (sum(I_active==i)==0) % i khong thuoc I active
    %         Am_X(i,i)=1;
    %         Am_Y(i,i)=1;
    %     end
    % end
    
    %tinh dp/dx, dp/dy
    pxy=Am\[-Am_X*vtp+Bm_X -Am_Y*vtp+Bm_Y];
    pxy_dot=Am\[Bm_X_dot Bm_Y_dot];
    
    pnx=pxy(:,1);
    pny=pxy(:,2);
    pnx_dot=pxy_dot(:,1);
    pny_dot=pxy_dot(:,2);
    
    D_jc=-[vt_S'*pnx vt_S'*pny;...
        vt_R'*pnx vt_R'*pny];
    
    
    D_jc_dot=-[vt_S'*pnx_dot vt_S'*pny_dot;...
        vt_R'*pnx_dot vt_R'*pny_dot];
    
    fx=-vt_S'*vtp
    fy=-vt_R'*vtp
    
    u_new=u-D_jc\([fx;fy]-[Wx;Wy])
    
    if ((norm(u_new-u)/norm(u_new))<=1e-7&&(norm([fx;fy]-[Wx;Wy])/norm([fx;fy]))<1e-7)
        disp('done')
        break
    end
    x_table(thu)=x;
    y_table(thu)=y;
    p_table(thu)=max(vtp)*6*nuy*omega*(R/c)^2;
    fx_table(thu)=fx*(6*nuy*omega*(R/c)^2*R*R);
    fy_table(thu)=fy*(6*nuy*omega*(R/c)^2*R*R);
    
    u=u_new;
    thu=thu+1
end
%%
table=[x_table' y_table' p_table' fx_table' fy_table'];
%%
disp('Pmax:'),max(vtp*6*nuy*omega*(R/c)^2+Pc)
cd=linspace(0,2,m+1)*L/2*1000;
mtp=zeros(m+1,n+1);
% mtq=zeros(m+1,n+1);
mtp(:,1:n)=vtp(luoi(:,1:n));
% mtq(:,1:n)=vtq(luoi(:,1:n));
mtp(:,n+1)=mtp(:,1);
% mtq(:,n+1)=mtq(:,1);
[X,Y]=meshgrid(cd,phi/pi*180);
%%
figure(1)
s=surf(X',Y',mtp);
axis tight
xlabel('Bearing length - L(mm)','fontweight','bold','fontsize',25,'color','black');
ylabel('Cicumference - \theta (rad)','fontweight','bold','fontsize',25,'color','black');
zlabel('Pressure - P(Pa)','fontweight','bold','fontsize',25,'color','black')
yticks([0:60:360])
xticks([0:10:50])
view(45,45)
alpha 0.8
pbaspect([1 1.5 1])
% i = rotate3d;set(i,'ActionPreCallback',...
%     'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(i,'ActionPostCallback',...
%     'set(gcf,''windowbuttonmotionfcn'','''')')
lighting phong

set(s,'edgecolor',[0 0 0.4],'meshstyle','both','linewidth',.15);
%%
figure(2)
t=mesh(X',Y',mtp*0)
xlabel('Bearing length - L(mm)','fontweight','bold','fontsize',25,'color','black');
ylabel('Cicumference - \theta (^o)','fontweight','bold','fontsize',25,'color','black');
yticks([0:60:360])
xticks([0:10:50])

%%
figure(3)
plot(phi/pi*180,h*c)
% axis([0 360 1.94e-5 2.06e-5])
xticks([0:60:360])
figure(4)
mtR(:,1:n)=vt_R(luoi(:,1:n));
mtS(:,1:n)=vt_S(luoi(:,1:n));
mtR(:,n+1)=mtR(:,1);
mtS(:,n+1)=mtS(:,1);
surf(X',Y',mtR)
figure(5)
surf(X',Y',mtS)
%%
figure(20)
pgiua=mtp(m/2+1,1:n+1)*6*nuy*omega*(R/c)^2/30e4;
r=5;
theta=(0:0.01:2*pi)+pi/2;
[x20,y20]=pol2cart(theta,r);
plot(x20,y20,'-k','LineWidth',2)
[x20,y20]=pol2cart(phi(1:n+1)+pi/2,r+pgiua);
maxcs=find(pgiua==max(pgiua));
[xmax,ymax]=pol2cart(phi(maxcs)+pi/2,r+pgiua(maxcs));

% line([xmax,0],[ymax,0])
hold on
% plot(xmax,ymax,'*')

plot(x20,y20,'-','LineWidth',2)
axis equal
line([u(2)*30,-u(2)*30],[-u(1)*30,u(1)*30])
theta=(0:0.01:2*pi)+pi/2;
r=4;
[x20,y20]=pol2cart(theta,r);
plot(x20+u(2)*3,y20-u(1)*3,'-','LineWidth',2)
plot(u(2)*3,-u(1)*3,'.','MarkerSize',10)
plot(xmax,ymax,'*')
line([0,0],[-6,6],'Color','black')
line([-6,6],[-0,0],'Color','black')
plot(0,0,'.k','MarkerSize',10)
line([u(2)*3,u(2)*3],[-6+-u(1)*3,6+-u(1)*3],'Color','black')
line([-6+u(2)*3,6+u(2)*3],[-0+-u(1)*3,0+-u(1)*3],'Color','black')
% axis([-8 8 -8 11])
r=5;
theta=(0:0.01:2*pi)+pi/2;
[x20,y20]=pol2cart(theta,r);
plot(x20,y20,'-k','LineWidth',2)
%%
figure(6)
o=1+(x*cos(phi)+y*sin(phi));
% yyaxis left
plot(phi(1:n+1)/pi*180,mtp(m/2+1,1:n+1)*6*nuy*omega*(R/c)^2,'k-');
hold on
xticks([0:60:360])
axis([0 360 0 9e4])
xlim([0,360])
plot(phi(1:n+1)/pi*180,mtp(11,1:n+1)*6*nuy*omega*(R/c)^2,'k-');
plot(phi(1:n+1)/pi*180,mtp(5,1:n+1)*6*nuy*omega*(R/c)^2,'k-');
hold off
xlabel('cicumference - \theta (^{o})','fontweight','bold','fontsize',11,'color','black');
ylabel('pressure - P(Pa)','fontweight','bold','fontsize',11,'color','black');

grid on 
%%
% yyaxis right
% plot(phi(1:n+1)/pi*180,o*c*1000,'k-');
% ylabel('film thickness  - h(mm)','fontweight','bold','fontsize',11,'color','black')
% ylim([0,0.07])


%%
% figure(7)
% plot(cd,mtp(:,bien_capdau))
% hold on
% plot(cd,mtp(:,bien_capdau+1))
% plot(cd,mtp(:,bien_capdau+2))
% plot(cd,mtp(:,bien_capdau+3))
% plot(cd,mtp(:,bien_capdau+4))
% hold off
% for i=1:length(bien_cd)
% [row,col]=find(bien_cd(i)==luoi);
% end
%%
figure(8)
Kxx=6*nuy*omega*R^2*(R/c)^2*D_jc(1,1)/c;
Kxy=6*nuy*omega*R^2*(R/c)^2*D_jc(1,2)/c;
Kyx=6*nuy*omega*R^2*(R/c)^2*D_jc(2,1)/c;
Kyy=6*nuy*omega*R^2*(R/c)^2*D_jc(2,2)/c;
K=[Kxx Kxy;Kyx Kyy];
% B{ll}=6*nuy*omega*R^2*(R/c)^2*D_jc_dot/(Wx*6*nuy*omega*(R/c)^2*R*R);
table=[x_table' y_table' p_table' fx_table' fy_table'];
% end

Kyx=abs(Kyx);
loglog([0.1:0.1:1 2:1:10],Kxx(:),'-',[0.1:0.1:1 2:1:10],Kxy(:),'-',[0.1:0.1:1 2:1:10],Kyx(:),'-',[0.1:0.1:1 2:1:10],Kyy(:),'-')
maxP=max(vtp*6*nuy*omega*(R/c)^2)