syms h_ct y_ct x_ct h_ct dh_dx_ct 
h_ct(x_ct,y_ct)=1 + (h1-h2)/h2*(1 - acot(y_ct/x_ct)/beta);%cong thuc h
dh_dx_ct(x_ct,y_ct)=diff(h_ct,x_ct);%dh/dx
