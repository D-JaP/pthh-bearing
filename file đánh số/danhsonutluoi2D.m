
function [luoi]=danhsonutluoi2D(kchia_m,kchia_n)% danh so nut
p_so=1;
luoi=zeros(kchia_m+1,kchia_n+1);


for j=1:kchia_n+1
    if mod(j,2)~=0
        for i=1:kchia_m+1
            luoi(i,j)=p_so;
            p_so=p_so+1;
        end
    end
    if mod(j,2)==0
        for i=1:kchia_m+1
            luoi(i,j)=p_so;
            p_so=p_so+1;
        end
    end
end