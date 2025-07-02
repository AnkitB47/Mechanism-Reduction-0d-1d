function obj = setTemp(obj,temp)
% Sets the temp vector, empty argument clears the 
% temps property
% obj = obj.setTemp(temps)
% obj = setTemp()
if nargin==2
    obj.temp = temp;
else
    obj.temp = [];
end
end