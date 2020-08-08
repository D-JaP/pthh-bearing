function luoi_phantu3D=luoiphantu3D(kchia_m,kchia_n,kchia_q)
luoi_phantu3D=zeros(kchia_m,kchia_n,kchia_q);
dem=1;
for k=1:kchia_q
    for j=1:kchia_n
        for i=1:kchia_m
            luoi_phantu3D(i,j,k)=dem;
            dem=dem+1;
        end
    end
end
