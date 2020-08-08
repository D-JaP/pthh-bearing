function J=Jctugiac(x,y,eta,xi)
J=1/4*[(1-eta)*(x(2)-x(1))+(1+eta)*(x(3)-x(4))   (1-eta)*(y(2)-y(1))+(1+eta)*(y(3)-y(4));...
       (1-xi )*(x(4)-x(1))+(1+xi )*(x(3)-x(2))   (1-xi )*(y(4)-y(1))+(1+xi )*(y(3)-y(2))];
end
