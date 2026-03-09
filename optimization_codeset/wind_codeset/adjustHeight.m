function [U_1] = adjustHeight(U_0,h_0,h_1,type,constant)

%see renewable energy resources page 294
%see wind energy explained equation 2.34
%updated 3/8/2026 to be able to run log type as a vector calc

if isequal(type,'power')
    U_1 = U_0*(h_1/h_0)^(constant);
elseif isequal(type,'log') 
    U_1 = U_0.*(log(h_1./(constant./1000))./log(h_0./(constant./1000)));
end

end

