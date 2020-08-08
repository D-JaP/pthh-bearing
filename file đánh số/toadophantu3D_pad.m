function Coordinate3D_polar_pad=toadophantu3D_pad(snut3D,kchia_m,kchia_n,kchia_q,r1,r2,beta,luoi_pad3D,khoang_i_pad,t_pad)
Coordinate3D_polar_pad=zeros(3,snut3D);
for i=khoang_i_pad
    [R,S,T]=ind2sub(size(luoi_pad3D),find(luoi_pad3D==i));
    Coordinate3D_polar_pad(1,i)=pi/2-(S-1)*beta/kchia_n;
    Coordinate3D_polar_pad(2,i)=(r2-(r2-r1)/kchia_m*(R-1))/r1;
    Coordinate3D_polar_pad(3,i)=(T-1)/kchia_q*t_pad/r1;
end