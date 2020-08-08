function t=chuyennutvaophantu(spt2D,luoi_phantu,kchia_m,kchia_n,luoi)
t=zeros(4,spt2D);
for i=1:kchia_m*kchia_n
    [row,col]=find(luoi_phantu==i);
%     if (col~=n)

        t(1:4,i)=[luoi((row+1),col),luoi((row+1),col+1),luoi(row,col+1),luoi(row,col)];
        
%     else
%         t(1:4,i)=[luoi((row+1),col),luoi((row+1),1),luoi(row,1),luoi(row,col)];
%     end
end