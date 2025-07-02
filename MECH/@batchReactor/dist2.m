function val = dist2(obj1,obj2,speciesList)
%val = dist2(obj1,obj2,speciesList)
% Simple distance function of two chemistry objects.
% It measures the relative differences in temperature  and in the species
% given by speciesList
% sqrt(sum(int(X_i(obj1)-X_i(obj2))^2))/int(X_i(obj1))
% Note that the norm is scaled by int(obj1). 

% speciesList = obj1.getInitalSpecies;

time = grefine(obj1.times);

n = length(speciesList);
if isa(speciesList,'cell')
    for k = 1:n
        indx(k) = getSpeciesIndex(speciesList{k});
    end
elseif isa(speciesList,'numeric')
     
    if min(speciesList)<=0
        error('dist2:ListContainsZero','speciesList should contain only nubers =>1')
    end
    indx = speciesList;
else 
    error('dist2:WronFormat','SpeciesList should be a cell or numeric vector')
end

for k = 1:n
    moleFrac1 = interp1(obj1.times,obj1.getMoleFractions(indx(k)),time);
    moleFrac2 = interp1(obj2.times,obj2.getMoleFractions(indx(k)),time);
    int(k) = 1/n*sqrt(trapez(time,(moleFrac1-moleFrac2).^2))/sqrt(trapez(time,moleFrac1.^2)); 
end

temp1 = interp1(obj1.times,obj1.temp,time);
temp2 = interp1(obj2.times,obj2.temp,time);

val = sqrt(trapez(time,(temp1-temp2).^2))/sqrt(trapez(time,temp1.^2)); 
for k = 1:n
    if int(k)<inf
        val = val+int(k);
    end
end
val = val+(obj2.nLivingSpecies/obj1.nLivingSpecies)^2; 

% local functions

    function grid = grefine(grid)
        % refines the discrete time
        grid = sort([grid 0.5*(grid(2:end)+grid(1:end-1))]);
    end

    function val = trapez(x,y)
        % trapezidual rule
         val = sum(0.5*(y(2:end)+y(1:end-1)).*(x(2:end)-x(1:end-1)));        
    end
end