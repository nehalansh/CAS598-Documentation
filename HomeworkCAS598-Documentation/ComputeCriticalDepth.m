function y_c = ComputeCriticalDepth(Q,section_type,a,b,theta,alpha,gravity)

y_c = zeros(size(Q));

for ind = 1:length(Q)
    Z_star = Q(ind)./(cos(theta))^.5 / (gravity/alpha)^.5;
    f_yc= @(y_c_var) get_section_property(y_c_var,section_type,a,b,'A').*...
        (get_section_property(y_c_var,section_type,a,b,'D')).^.5 - Z_star;
    y_c_init = 3; 
    y_c(ind) = fzero(f_yc, y_c_init);
end

