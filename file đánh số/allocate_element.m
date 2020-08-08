function [I_active,I_boundry]=allocate_element(kchia_m,kchia_n,luoi)
dem_ib=1;
dem_ia=1;
for i=1:kchia_m+1
    for j=1:kchia_n+1
        
        if (i==1||j==1||i==kchia_m+1||j==kchia_n+1)
            I_boundry(dem_ib)=luoi(i,j);
            dem_ib=dem_ib+1;
        else
            I_active(dem_ia)=luoi(i,j);
            dem_ia=dem_ia+1;
        end
    end
end