function obj = setDens(obj,dens)
% Sets the dens vector, empty argument clears the 
% times property
% obj = obj.setDens(times)
% obj = setDens()
if nargin==2
    obj.dens = dens;
else
    obj.dens = [];
end
end