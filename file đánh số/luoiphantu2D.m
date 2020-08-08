function luoi_phantu=luoiphantu2D(kchia_m,kchia_n)
p_so2=1;
luoi_phantu=zeros(kchia_m,kchia_n);
for j=1:kchia_n
    %     if mod(i,2)~=0
    for i=1:kchia_m
        luoi_phantu(i,j)=p_so2;
        p_so2=p_so2+1;
    end
    %     end
    %     if mod(i,2)==0
    %         for j=2*n:-1:1
    %             luoi_phantu(j,i)=p_so2;
    %             p_so2=p_so2+1;
    %         end
    %     end
end