function y_c = ComputeCriticalDepth(Q,section_type,a,b,theta,alpha,gravity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% y_c = ComputeCriticalDepth(Q,section_type,a,b,theta,alpha,gravity)
%
% This function computes the critical depth, y_c, for a certain flow in channel
% sections of different shapes.
%
% INPUTS:
%
% Q = flow. It can be an array of values.
%
% section_type = string that can be either: 'trapezoidal', 'rectangular',
% 'triangular', 'generic':
%
% a, b = paramaters desribing the section. See the help of the function
% "get_section_property.m".
%
% theta = angle of slope of the channel bottom;
%
% alpha: kinetic energy correction factor
%
% gravity = value of the gravity acceleration
%
% OUTPUT:
%
% y_c = critical depth correspoinding to Q (it can be an array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_c = zeros(size(Q));

for ind = 1:length(Q)

    % We need to solve this equation:
    % A*(y_c)^.5 = Q/(cos(theta))^.5 / (g/alpha)^.5
    
    % Compute constant
    Z_star = Q(ind)./(cos(theta))^.5 / (gravity/alpha)^.5;
    
    % Define function as equation
    f_yc= @(y_c_var) get_section_property(y_c_var,section_type,a,b,'A').*...
        (get_section_property(y_c_var,section_type,a,b,'D')).^.5 - Z_star;

    % Solve equation
    y_c_init = 3; 
    y_c(ind) = fzero(f_yc, y_c_init);
end

