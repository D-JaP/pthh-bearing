function [xi,eta,W]=tpso_Gauss_2d(a,demdiemtichphan)
switch a
    case 1
    case 2
        switch demdiemtichphan
            case 1
                xi=-1/sqrt(3);
                eta=-1/sqrt(3);
                W(demdiemtichphan)=1;
            case 2
                xi=1/sqrt(3);
                eta=-1/sqrt(3);
                W(demdiemtichphan)=1;
            case 3
                xi=1/sqrt(3);
                eta=1/sqrt(3);
                W(demdiemtichphan)=1;
            case 4
                xi=-1/sqrt(3);
                eta=1/sqrt(3);
                W(demdiemtichphan)=1;    
        end
    case 3
        switch demdiemtichphan
            case 1 
                xi=-sqrt(0.6);
                eta=-sqrt(0.6);
                rutgontrongluong3diemtichphan
            case 2
                xi=0;
                eta=-sqrt(0.6);
                rutgontrongluong3diemtichphan
            case 3
                xi=sqrt(0.6);
                eta=-sqrt(0.6);
                rutgontrongluong3diemtichphan
            case 4
                xi=-sqrt(0.6);
                eta=0;
                rutgontrongluong3diemtichphan
            case 5
                xi=0;
                eta=0;
                rutgontrongluong3diemtichphan
            case 6
                xi=sqrt(0.6);
                eta=0;
                rutgontrongluong3diemtichphan
            case 7
                xi=-sqrt(0.6);
                eta=sqrt(0.6);
                rutgontrongluong3diemtichphan
            case 8
                xi=0;
                eta=sqrt(0.6);
                rutgontrongluong3diemtichphan
            case 9
                xi=sqrt(0.6);
                eta=sqrt(0.6);
                rutgontrongluong3diemtichphan
        end
end