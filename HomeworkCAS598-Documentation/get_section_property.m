 function Output_VAR = get_section_property(y,section_type,a,b,property)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output_VAR = get_section_property(y,section_type,a,b,property)
%
% This script calculates one of these outputs:
%
% area (A), hydraulic radius (R), wet perimeter (P),
% top width (B), and hydraulic depth (D) for different types of section,
%
% given the flow depth y and the section properties.
% The type of output is specified by
% property = string that can be either: 'A', 'R', 'P', 'B', 'D' depending on
% the desired property.
%
% INPUTS:
%
% section_type = string that can be either:
%
% 'trapezoidal':
% a = channel width at the bottom
% b = horizontal length for 1 vertical rise
%
%                           \          /
%                            \        /   |1
%                             \------/ -b-
%                             <--a-->
%
% 'rectangular':
% a = bottom width
% b = []
%
% 'circular':
% a = diameter
% b = []
%
% 'triangular':
% a = horizontal length for 1 vertical rise
%
%                           \     /
%                            \   /   | 1
%                             \ / -a-
%                              V
%
% 'generic':
% a = array with x coordinates of the points surveyed in the cross-section
% b = array with y coordinates of the points surveyed in the cross-section
% NOTE: points have to be inserted from left to right adopting the
% convention that the left (right) bank is located to the left (right) of
% an observer that is looking downstream.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch section_type
    
    case 'trapezoidal'
        A = (a + b.*y).*y;
        P = (a + 2.*y.*(1+b^2)^.5);
        B = a + 2*b.*y;
        
    case 'rectangular'
        A = a .* y;
        P = a + 2.*y;
        B = a;
        
    case 'circular'
        D = a;
        Fi = 2*acos(1 - y./(D/2));
        A = (D/2)^2/2*(Fi - sin(Fi));
        P = D/2*Fi;
        B = D*sin(Fi/2);
        
    case 'triangular'
        A = y.^2.*a;
        P = 2.*y.*(1 + a.^2).^.5;
        B = 2.*a.*y;
        
    case 'generic'
        
        x_section = a; y_section = b;
        %         % If there are vertical parts, add a small number to eliminate them
        %         dum = find(diff(x_section) == 0);
        %         if ~isempty(dum)
        %             for ind = 1:length(dum)
        %                 x_section(dum(ind)) = x_section(dum(ind)) - 10^(-6);
        %             end
        %         end
        
        % Initialize three vectors with same size as y
        A = NaN*ones(size(y));
        P = NaN*ones(size(y));
        B = NaN*ones(size(y));
        
        for ind_y = 1:length(y)
            
            y_star = y(ind_y) - y_section;
            
            % Identify the points that we want to keep
            ind_to_keep = find(y_star >= 0);
            x_integr = x_section(ind_to_keep);
            y_integr = y_section(ind_to_keep);
            
            % If needed, add points through interpolation
            % Find the indices used for the interpolation
            dum = y_star >= 0; dum = diff(dum);
            % Find the indices of x_section and y_section that are
            % immediately above the water depth y
            IND_LEFT = find(dum == 1);
            IND_RIGHT = find(dum == -1)+1;
            
            % Note that we can have more than one pair where we need to
            % interpolate beacue the section can be above the level of the
            % water in the middle of the river (e.g., case of an island)
            % Initialize a matrix with two columns (one for the left and right
            % x coordinates) and number of rows equal to the number of
            % pairs
            x_interp = NaN(length(IND_LEFT),2);
            y_interp = ones(length(IND_LEFT),2)*y(ind_y);
            for ind_set = 1:length(IND_LEFT)
                ind_dum1 = IND_LEFT(ind_set);
                ind_dum2 = ind_dum1+1;
                x_interp(ind_set,1) = interp1(y_section([ind_dum1 ind_dum2]),x_section([ind_dum1 ind_dum2]),y(ind_y),'linear');
                ind_dum1 = IND_RIGHT(ind_set);
                ind_dum2 = ind_dum1-1;
                x_interp(ind_set,2) = interp1(y_section([ind_dum1 ind_dum2]),x_section([ind_dum1 ind_dum2]),y(ind_y),'linear');
            end
            
            % If we have more than pair, focus only on the one with the
            % largest top width
            
            % Compute B as the largest width
            [B, IND_MAX] = max(diff(x_interp,[],2));
            
            % Include the interpolated values in x_integr and y_integr
            x_integr = [x_interp(IND_MAX,:) x_integr];
            y_integr = [y_interp(IND_MAX,:) y_integr];
            
            [x_integr, ind_sorted] = sort(x_integr);
            y_integr = y_integr(ind_sorted);
            
            % Compute area using the integral
            y_star = y(ind_y) - y_integr;
            A(ind_y) = trapz(x_integr,y_star);
            
            % Get the wet perimeter
            % For this aim, let us get the differences of x_integr and
            % y_integr. These are the length of triangle edges.
            % Then, we need to compute the hypothenuse.
            diff_x_dum = diff(x_integr);
            diff_y_dum = diff(y_integr);
            hypotenuse = (diff_x_dum.^2 + diff_y_dum.^2).^.5;
            
            % Nel vettore diff_y_dum, i valori pari a zero non contribuiscono al
            % contorno bagnato ---------  Why??????
            hypotenuse(diff_y_dum == 0 & y_star(2:end) == 0) = 0;
            
            P(ind_y) = sum(hypotenuse);
            
            figure
            plot(x_section,y_section,'ko-'), hold on
            plot([min(x_section) max(x_section)],[y(ind_y) y(ind_y)],'--','Linewidth',2)
            plot(x_integr,y_integr,'rx')
        end
        
    otherwise
        error('Wrong section type! It can be on of these: ''trapezoidal'', ''rectangular'', ''circular'', ''triangular'', or ''generic''');
        
end

R = A./P;
D = A./B;

% Assign the output
switch property
    case 'A'
        Output_VAR = A;
    case 'R'
        Output_VAR = R;
    case 'B'
        Output_VAR = B;
    case 'P'
        Output_VAR = P;
    case 'D'
        Output_VAR = D;
    otherwise
        error('The property %s is not available!',property);
end



