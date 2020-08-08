if (eta~=0&&xi~=0)
    W(demdiemtichphan)=5/9*5/9;   
elseif ((eta~=0&&xi==0)||(eta==0&&xi~=0))
    W(demdiemtichphan)=5/9*8/9;
elseif (eta==0&&xi==0)
    W(demdiemtichphan)=8/9*8/9;
end