function d_dr_bar=d_dr_bar(vecto_candaoham)
global kchia_m kchia_n kchia_q luoi3D r1 r2 snut3D
d_dr_bar=zeros(snut3D,1);
for k=1:kchia_q+1
    for j=1:kchia_n+1
        for i=kchia_m+1:-1:2
            d_dr_bar(luoi3D(i,j,k))=(vecto_candaoham(luoi3D(i-1,j,k))-vecto_candaoham(luoi3D(i,j,k)))/(r2-r1)*kchia_m*r1;
        end
        d_dr_bar(luoi3D(1,j,k))=interp1(2:kchia_m+1,d_dr_bar(luoi3D(2:kchia_m+1,j,k)),1,'linear','extrap');
    end
end