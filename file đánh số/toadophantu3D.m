function Coordinate3D_polar=toadophantu3D(snut3D,kchia_m,kchia_n,kchia_q,r1,r2,beta,H3D,H3D_global,luoi3D)
Coordinate3D_polar=zeros(3,snut3D);
for i=1:snut3D
    [R,S,T]=ind2sub(size(luoi3D),find(luoi3D==i));
    Coordinate3D_polar(1,i)=(S-1)/kchia_n;
    Coordinate3D_polar(3,i)=(r2-(r2-r1)/kchia_m*(R-1))/r1;
    Coordinate3D_polar(2,i)=H3D(i)/H3D_global(i);
end