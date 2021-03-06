function [mtA,f]=tichphanbien2D(index_phantu, so_phan_tu,toado_phantu, T_alpha, snut_all,check_1_2)
global r1 r2 Co_h Co_k

mtA=zeros(snut_all);
f=zeros(snut_all,1);
    for i=1:so_phan_tu 
            xP=toado_phantu(1,index_phantu(1:4,i));
            yP=toado_phantu(2,index_phantu(1:4,i));
        for j=1:4
            [xi,eta,trong_so]=tpso_Gauss_2d(2,j);
            [N1,N2,N3,N4]=hamdang(xi,eta);
            N=[N1 ;N2 ;N3 ;N4];
            dNij=[ -1/4*(1-eta) -1/4*(1-xi);...
                1/4*(1-eta) -1/4*(1+xi);...
                1/4*(1+eta)  1/4*(1+xi);...
                -1/4*(1+eta)  1/4*(1-xi)];
            
            J=Jctugiac(xP,yP,eta,xi);
            Jn=inv(J);
            
%             dNx(1:4)=Jn(1,1)*dNij(1:4,1)+Jn(1,2)*dNij(1:4,2);
%             dNz(1:4)=Jn(2,1)*dNij(1:4,1)+Jn(2,2)*dNij(1:4,2);
            
            
            A=N'*xP';%toa do de cac
            B=N'*yP';%toa do de cac
            
            %tinh cac ma tran do cung
            switch check_1_2
                case 1 %r1
                    A_i=N*N'.*(Co_h*r1/Co_k).*(r2/r1).*abs(det(J)).*trong_so(j);%toa do cuc
                    mtA(index_phantu(1:4,i),index_phantu(1:4,i))=mtA(index_phantu(1:4,i),index_phantu(1:4,i))+A_i;

                    f_i=N'*(Co_h*r1/Co_k)*(r2/r1)*T_alpha*abs(det(J))*trong_so(j);%toa do cuc
                    f(index_phantu(1:4,i),1)=f(index_phantu(1:4,i),1)+f_i';
                case 2  %r2
                    A_i=N*N'.*(Co_h*r1/Co_k).*(r1/r1).*abs(det(J)).*trong_so(j);%toa do cuc
                    mtA(index_phantu(1:4,i),index_phantu(1:4,i))=mtA(index_phantu(1:4,i),index_phantu(1:4,i))+A_i;

                    f_i=N'*(Co_h*r1/Co_k)*(r1/r1)*T_alpha*abs(det(J))*trong_so(j);%toa do cuc
                    f(index_phantu(1:4,i),1)=f(index_phantu(1:4,i),1)+f_i';
                case 3  %z1
                    A_i=N*N'.*(Co_h*r1/Co_k).*(B/r1).*abs(det(J)).*trong_so(j);%toa do cuc
                    mtA(index_phantu(1:4,i),index_phantu(1:4,i))=mtA(index_phantu(1:4,i),index_phantu(1:4,i))+A_i;
                    
                    f_i=N'*(Co_h*r1/Co_k)*(B/r1)*T_alpha*abs(det(J))*trong_so(j);%toa do cuc
                    f(index_phantu(1:4,i),1)=f(index_phantu(1:4,i),1)+f_i';
                case 4 
                    A_i=N*N'.*(Co_h*r1/Co_k)/(A/r1).*abs(det(J)).*trong_so(j);%toa do cuc
                    mtA(index_phantu(1:4,i),index_phantu(1:4,i))=mtA(index_phantu(1:4,i),index_phantu(1:4,i))+A_i;
                    
                    f_i=N'*(Co_h*r1/Co_k)/(A/r1)*T_alpha*abs(det(J))*trong_so(j);%toa do cuc
                    f(index_phantu(1:4,i),1)=f(index_phantu(1:4,i),1)+f_i';
                case 5
                    A_i=N*N'.*(Co_h*r1/Co_k)/(A/r1).*abs(det(J)).*trong_so(j);%toa do cuc
                    mtA(index_phantu(1:4,i),index_phantu(1:4,i))=mtA(index_phantu(1:4,i),index_phantu(1:4,i))+A_i;
                    
                    f_i=N'*(Co_h*r1/Co_k)/(A/r1)*T_alpha*abs(det(J))*trong_so(j);%toa do cuc
                    f(index_phantu(1:4,i),1)=f(index_phantu(1:4,i),1)+f_i';
            end
        end
    end
    