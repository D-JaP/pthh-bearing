function [luoi3D]=danhsonutluoi3D(kchia_m,kchia_n,kchia_q)% danh so nut
p_so=1; % bien chay
luoi3D=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
for k=1:kchia_q+1
    for j=1:kchia_n+1
        if mod(j,2)~=0
            for i=1:kchia_m+1
                luoi3D(i,j,k)=p_so;
                p_so=p_so+1;
            end
        end
        if mod(j,2)==0
            for i=1:kchia_m+1
                luoi3D(i,j,k)=p_so;
                p_so=p_so+1;
            end
        end
    end
end