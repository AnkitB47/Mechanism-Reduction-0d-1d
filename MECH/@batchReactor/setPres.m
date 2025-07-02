function obj = setPres(obj,pres)
% Sets the pres vector, empty argument clears the 
% pres property
% obj = obj.setPres(times)
% obj = setPres()
if nargin==2
    obj.pres = pres;
else
    obj.pres = [];
end
end