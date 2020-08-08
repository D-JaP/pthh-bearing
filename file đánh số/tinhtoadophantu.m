function p=tinhtoadophantu(snut,luoi,kchia_m,kchia_n,beta,r1,r2)
p=zeros(2,snut);
for i=1:snut
    [row,col]=find(luoi==i);
    p(1,i)=pi/2-(col-1)*beta/kchia_n;
    p(2,i)=(r2-(r2-r1)/kchia_m*(row-1));
end
