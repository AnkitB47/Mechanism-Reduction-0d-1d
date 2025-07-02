function CanteraString = vector2CanteraString(speciesNames,v)
% CanteraString = vector2CanteraString(speciesNames,vectorOfFractions)
% Function to create a Mass/Molefractions string for use with canteras
% set(gas) method. Useful to initialize a gas object with values computed by
% onother program or given data.
% Exapmle:
% Assume we have a chemistry object c with a solved reactor
% We want to use send the mole fraction at the end of the time domain to a
% gas object.
%   moles = c.getMoleFractions({'CH4' 'O2'});
%   set(gas,'X',vector2CanteraString({'CH4' 'O2'},moles(:,end)));
%  (c) Hong Bin Xu for VIRTUHCON


if ~iscellstr(speciesNames)
    error('Vector2CanteraString:MustBeCell','The Names of the species must be a cell array of strings');
end
if ~isnumeric(v)
    error('Vector2CanteraString:NotNumeric','The array of the fractions must be numeric');
else
    if length(speciesNames)~=length(v)
        error('Vector2CanteraString:DimensionMissmatch','The dimensions of the names and the fractions must be equal');
    end
end

CanteraString = '';
n = length(speciesNames);
for i = 1:n
    S = char(speciesNames(i));
    CanteraString = strcat(CanteraString,S,':',num2str(v(i)));
    if i~=n
        CanteraString=strcat(CanteraString,',');
    end
end
end