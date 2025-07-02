function equilibrate(obj,arg)
% varargout = obj.equilibrate  
% Additional argument must be a two-letter string, which must
% be one of the set
%        ['TP','TV','HP','SP','SV','UV','PT','VT','PH','PS','VS','VU'].
%        If H, U, S, or V is specified, the value must be the specific
%        value (per unit mass).
 
if isempty(obj.initialMoleFractions)&&isempty(obj.initialMassFractions)
    throw(emptyFieldException)
end
if isempty(obj.initialMoleFractions)
    mole = false;
else
    mole = true;
end
fprintf('Compute mole fractions y/n: ');
if mole
    fprintf('y\n')
else
    fprintf('n\n')
end
if nargin<2
    arg = 'TP';
end
try
    set(obj.gas,'T',obj.initialTemperature);
    set(obj.gas,'P',obj.initialPressure);
    if mole
        set(obj.gas,'X',obj.initialMoleFractions);
    else
        set(obj.gas,'Y',obj.initialMassFractions);
    end
     
    equilibrate(obj.gas,arg);
     
catch ME
    throw(ME)
end

switch nargout
    case 0
        obj.temp = temperature(obj.gas);
        obj.pres = pressure(obj.gas);
        obj.dens = density(obj.gas);
        if mole
            obj.moleFractions = moleFractions(obj.gas);
        else
            obj.massFractions = massFractions(obj.gas);
        end
    otherwise
        throw(wrongNumberOutputException);
end
end     
