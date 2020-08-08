function [xi,W]=tpso_Gauss_1d(a,demdiemtichphan)
switch a
    case 1
    case 2
        switch demdiemtichphan
            case 1
                xi=-1/sqrt(3);
                W(demdiemtichphan)=1;
            case 2
                xi=1/sqrt(3);   
                W(demdiemtichphan)=1;
        end
    case 3
        switch demdiemtichphan
            case 1 
                xi=-sqrt(0.6);
                W(demdiemtichphan)=5/9;
            case 2
                xi=0;
                W(demdiemtichphan)=8/9;
            case 3
                xi=sqrt(0.6);
                W(demdiemtichphan)=5/9;
        end
end