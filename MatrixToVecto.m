function [vecto]=MatrixToVecto(M)
[a,b]=size(M);
vecto=zeros(a*b,1);
q=1;
for i=1:a
    for j=1:b
        vecto(q)=M(i,j);
        q=q+1;
    end
end
