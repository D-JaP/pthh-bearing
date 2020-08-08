function [mtA,mtAh,dA_dphi,dA_dr,db_dphi,db_dr,b,bh,unit_vt,vt_X,vt_Y]=tichphanGauss2D_ochan_nhiet()
global snut spt2D  t pp p delta beta h1 h2 nuy omega
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