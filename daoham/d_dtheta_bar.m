function d_dtheta_bar=d_dtheta_bar(vecto_candaoham)
global kchia_m kchia_n kchia_q luoi3D beta snut3D
d_dtheta_bar=zeros(snut3D,1);
for k=1:kchia_q+1
    for i=1:kchia_m+1
        for j=2:kchia_n+1
            d_dtheta_bar(luoi3D(i,j,k))=(vecto_candaoham(luoi3D(i,j,k))-vecto_candaoham(luoi3D(i,j-1,k)))/beta*kchia_n;
        end
        d_dtheta_bar(luoi3D(i,1,k))=interp1(2:kchia_n+1,d_dtheta_bar(luoi3D(i,2:kchia_n+1,k)),1,'linear','extrap');
    end
end