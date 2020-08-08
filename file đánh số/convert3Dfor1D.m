function [convert3Dfor1D]=convert3Dfor1D(kchia_m,kchia_n,kchia_q,luoi3D,a)% danh so nut
convert3Dfor1D=zeros((kchia_m+1)*(kchia_n+1)*(kchia_q+1),(kchia_q+1));
for k=1:kchia_q+1
    for j=1:kchia_n+1
        for i=1:kchia_m+1
            convert3Dfor1D(luoi3D(i,j,k),:)=a(luoi3D(i,j,1:kchia_q+1));
        end
    end
end