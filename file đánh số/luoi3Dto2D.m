function [luoi3Dto2D]=luoi3Dto2D(kchia_m,kchia_n,kchia_q,luoi)% danh so nut
luoi3Dto2D=zeros(kchia_m+1,kchia_n+1,kchia_q+1);
for k=1:kchia_q+1
    for j=1:kchia_n+1
        for i=1:kchia_m+1
            luoi3Dto2D(i,j,k)=luoi(i,j);
        end
    end
end