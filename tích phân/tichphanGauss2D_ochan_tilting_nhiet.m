tichphanGauss2D_ochan_tilting_nhiet
% global snut spt2D t pp p nuy omega alpha_r alpha_theta hp
    mtA=zeros(snut,snut);

%     mt_dA_dhp=zeros(snut,snut);% Partial derivatives A/hp
    mt_dA_dalpha_r=zeros(snut,snut);% Partial derivatives A/alpha_r
    mt_dA_dalpha_theta=zeros(snut,snut);% Partial derivatives A/alpha_theta
    
    b=zeros(snut,1);
%     bh=zeros(snut,1);
%     dB_dhp=zeros(snut,1);% Partial derivatives B/hp
    dB_dalpha_r=zeros(snut,1);% Partial derivatives B/alpha_r
    dB_dalpha_theta=zeros(snut,1);% Partial derivatives B/alpha_theta
    
    unit_vt=zeros(snut,1);%unit vecto for Gauss integral
    
    vt_X=zeros(snut,1);%unit vecto X for Gauss integral
    vt_Y=zeros(snut,1);%unit vecto Y for Gauss integral
    for i=1:spt2D
            xP=pp(1,t(1:4,i));
            yP=pp(2,t(1:4,i));
            theta=p(1,t(1:4,i));
            r=p(2,t(1:4,i));
        for j=1:4
            [xi,eta,W]=tpso_Gauss_2d(2,j);
            [N1,N2,N3,N4]=hamdang(xi,eta);
            N=[N1 ;N2 ;N3 ;N4];
            dNij=dN_2D(eta,xi);% tinh dNij
            
            J=Jctugiac(xP,yP,eta,xi);
            Jn=inv(J);
            
            dNx(1:4)=Jn(1,1)*dNij(1:4,1)+Jn(1,2)*dNij(1:4,2);
            dNz(1:4)=Jn(2,1)*dNij(1:4,1)+Jn(2,2)*dNij(1:4,2);
            
            A=N'*xP';
            B=N'*yP';
            C=N'*theta';%bien
            D=N'*r';%bien
            
            HH=N'*(H(t(1:4,i)))';
            HH2=HH^2;
            HH3=HH^3;
            
            dH_dalpha_theta=    eval(dh_dalpha_theta_ct(xP,yP,alpha_r, alpha_theta, hp));
            dH_dalpha_r=        eval(dh_dalpha_r_ct(xP,yP,alpha_r, alpha_theta, hp));
            %tinh cac ma tran do cung
            A_i=HH3                                          *(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
            dA_dalpha_r_i=3*HH2*dH_dalpha_r                  *(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
            dA_dalpha_theta_i=3*HH2*dH_dalpha_theta          *(dNx'*dNx+(dNz'*dNz))*abs(det(J))*W(j);
            
            mtA(t(1:4,i),t(1:4,i))=                            mtA(t(1:4,i),t(1:4,i))           	+A_i;
            mt_dA_dalpha_theta(t(1:4,i),t(1:4,i))=             mt_dA_dalpha_theta(t(1:4,i),t(1:4,i))  +dA_dalpha_theta_i;
            mt_dA_dalpha_r(t(1:4,i),t(1:4,i))=                 mt_dA_dalpha_r(t(1:4,i),t(1:4,i))      +dA_dalpha_r_i;
            %tinh vecto tai
%             A=xP;
%             B=yP;
%             dH=h2*delta*(B./(A.^2.*(B.^2./A.^2 + 1).*sqrt(A.^2 + B.^2)) + acot(B./A).*A./(A.^2 + B.^2).^(3/2));
            
            
            r=sqrt(A^2+B^2);
            
            dh_dx=eval(dh_dx_ct(xP,yP,alpha_r, alpha_theta, hp));
            dh_dy=eval(dh_dy_ct(xP,yP,alpha_r, alpha_theta, hp));
            Bj=-(6*nuy*omega*N*dh_dx*abs(det(J))*W(j)*sin(C)...
               -6*nuy*omega*N*dh_dy*abs(det(J))*W(j)*cos(C));
            
            %bhi=-6*nuy*omega*N*(A/sqrt(A^2 + B^2))*abs(det(J))*W(j)/h2^2;
%             bhi=-6*nuy*omega*N*B*(3*h1 - 2*h2)/(sqrt(A^2 + B^2)*h2^4*beta)*abs(det(J))*W(j);
            dB_dalpha_r_ii=eval(dB_dalpha_r_ct(xP,yP,alpha_r, alpha_theta, hp,C));
            dB_dalpha_theta_ii=eval(dB_dalpha_theta_ct(xP,yP,alpha_r, alpha_theta, hp,C));
            
            dB_dalpha_r_i=-6*nuy*omega*N*dB_dalpha_r_ii*abs(det(J))*W(j);
            dB_dalpha_theta_i=-6*nuy*omega*N*dB_dalpha_theta_ii*abs(det(J))*W(j);
            
            b(t(1:4,i))=b(t(1:4,i))+Bj;
%             bh(t(1:4,i))=bh(t(1:4,i))+bhi;
            dB_dalpha_r(t(1:4,i))=dB_dalpha_r(t(1:4,i))+dB_dalpha_r_i;
            dB_dalpha_theta(t(1:4,i))=dB_dalpha_theta(t(1:4,i))+dB_dalpha_theta_i;
            
            unit_vti=N*abs(det(J))*W(j);
            unit_vt(t(1:4,i))=unit_vt(t(1:4,i)) + unit_vti;
            
            X=N*A*abs(det(J))*W(j);
            Y=N*B*abs(det(J))*W(j);
            
            vt_X(t(1:4,i))=vt_X(t(1:4,i)) + X;%vecto x cho chuyen doi nodal forces
            vt_Y(t(1:4,i))=vt_Y(t(1:4,i)) + Y;%vecto y cho chuyen doi nodal forces
            
        end
    end