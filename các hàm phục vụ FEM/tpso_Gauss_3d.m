function [r,s,t,W]=tpso_Gauss_3d(a,demdiemtichphan)
switch a
    case 1
    case 2
        switch demdiemtichphan
            case 1
                r=-1/sqrt(3);
                s=-1/sqrt(3);
                t=-1/sqrt(3);
                W(demdiemtichphan)=1;
            case 2
                r=1/sqrt(3);
                s=-1/sqrt(3);
                t=-1/sqrt(3);
                W(demdiemtichphan)=1;
            case 3
                r=1/sqrt(3);
                s=1/sqrt(3);
                t=-1/sqrt(3);
                W(demdiemtichphan)=1;
            case 4
                r=-1/sqrt(3);
                s=1/sqrt(3);
                t=-1/sqrt(3);
                W(demdiemtichphan)=1;    
            case 5
                r=-1/sqrt(3);
                s=-1/sqrt(3);
                t=1/sqrt(3);
                W(demdiemtichphan)=1;
            case 6
                r=1/sqrt(3);
                s=-1/sqrt(3);
                t=1/sqrt(3);
                W(demdiemtichphan)=1;
            case 7
                r=1/sqrt(3);
                s=1/sqrt(3);
                t=1/sqrt(3);
                W(demdiemtichphan)=1;
            case 8
                r=-1/sqrt(3);
                s=1/sqrt(3);
                t=1/sqrt(3);
                W(demdiemtichphan)=1;    
        end
%     case 3
%         switch demdiemtichphan
%             case 1 
%                 r=-sqrt(0.6);
%                 s=-sqrt(0.6);
%                 t=-sqrt(0.6);
%                 rutgontrongluong3diemtichphan
%             case 2
%                 xi=0;
%                 eta=-sqrt(0.6);
%                 rutgontrongluong3diemtichphan
%             case 3
%                 xi=sqrt(0.6);
%                 eta=-sqrt(0.6);
%                 rutgontrongluong3diemtichphan
%             case 4
%                 xi=-sqrt(0.6);
%                 eta=0;
%                 rutgontrongluong3diemtichphan
%             case 5
%                 xi=0;
%                 eta=0;
%                 rutgontrongluong3diemtichphan
%             case 6
%                 xi=sqrt(0.6);
%                 eta=0;
%                 rutgontrongluong3diemtichphan
%             case 7
%                 xi=-sqrt(0.6);
%                 eta=sqrt(0.6);
%                 rutgontrongluong3diemtichphan
%             case 8
%                 xi=0;
%                 eta=sqrt(0.6);
%                 rutgontrongluong3diemtichphan
%             case 9
%                 xi=sqrt(0.6);
%                 eta=sqrt(0.6);
%                 rutgontrongluong3diemtichphan
%         end
end