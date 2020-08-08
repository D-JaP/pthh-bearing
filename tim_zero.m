function tim_zero(a)
dem=0;
[m,n]=size(a);
for i=1:m
    if (sum(a(i,:))==0)
        for j=1:n
            if (sum(a(:,j))==0)
%                 disp('co hang cot zero khien cho ma tran singular');
                dem=dem+1;
%                 disp(i)
%                 disp(j)
            end
        end
    end
end
dem