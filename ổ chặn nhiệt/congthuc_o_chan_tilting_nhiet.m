syms h_ct y_ct x_ct h_ct dh_dx_ct dh_dy_ct dh_dalpha_theta_ct dh_dalpha_r_ct alpha_r_ct alpha_theta_ct hp_ct dB_dalpha_r_ct dB_dalpha_theta_ct C
assume(x_ct,'positive')
assume(y_ct,'positive')

h_ct(x_ct,y_ct,alpha_r_ct, alpha_theta_ct, hp_ct) = 1 + sqrt(x_ct^2 + y_ct^2)*sin(theta_p - atan(x_ct/y_ct))*sin(alpha_r_ct)/hp_ct + (r_p - sqrt(x_ct^2 + y_ct^2)*cos(theta_p - atan(x_ct/y_ct)))*sin(alpha_theta_ct)/hp_ct;

dh_dx_ct(x_ct,y_ct,alpha_r_ct, alpha_theta_ct, hp_ct)=(diff(h_ct,x_ct));%dh/dx
dh_dy_ct(x_ct,y_ct,alpha_r_ct, alpha_theta_ct, hp_ct)=(diff(h_ct,y_ct));%dh/dx

dh_dalpha_theta_ct(x_ct,y_ct,alpha_r_ct, alpha_theta_ct, hp_ct)=diff(h_ct,alpha_theta_ct);
dh_dalpha_r_ct(x_ct,y_ct,alpha_r_ct, alpha_theta_ct, hp_ct)=diff(h_ct,alpha_r_ct);

dB_dalpha_r_ct(x_ct,y_ct,alpha_r_ct, alpha_theta_ct, hp_ct,C)=diff(dh_dx_ct*sin(C)-dh_dy_ct*cos(C),alpha_r_ct);
dB_dalpha_theta_ct(x_ct,y_ct,alpha_r_ct, alpha_theta_ct, hp_ct,C)=diff(dh_dx_ct*sin(C)-dh_dy_ct*cos(C),alpha_theta_ct);