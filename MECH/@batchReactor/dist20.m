function val = dist20(obj1,obj2)
speciesList = obj1.mainSpecies;
nsp = length(speciesList);
fprintf('\nCompute error for species:\n')
class( obj1.speciesNames(speciesList(1)))
for k = 1:nsp-1
    fprintf([ obj1.speciesNames{speciesList(k)},', ']);
    
end
fprintf([obj1.speciesNames{speciesList(k)},'\n']);

time = grefine(obj1.times);
time = grefine(time);

for k = 1:nsp            
    moleFrac1 = interp1(obj1.times,...
        obj1.massFractions(speciesList(k),:),time);
    moleFrac2 = interp1(obj2.times,...
        obj2.massFractions(speciesList(k),:),time);
    int(k) = 1/nsp*sqrt(trapez(time,(moleFrac1-moleFrac2).^2))...
        /sqrt(trapez(time,moleFrac1.^2)); 
end
temp1 = interp1(obj1.times,obj1.temp,time);
temp2 = interp1(obj2.times,obj2.temp,time);
val = sum(int) + ...
    trapez(time,(temp1-temp2).^2)/trapez(time,temp1.^2);
fprintf(['Relative error in species = ',num2str(val),'\n']);

if val>1
    val = Inf;
end
end 

% local functions

    function grid = grefine(grid)
        % refines the discrete time
        grid = unique([grid 0.5*(grid(2:end)+grid(1:end-1))]);
    end

    function val = trapez(x,y)
        % trapezidual rule
         val = sum(0.5*(y(2:end)+y(1:end-1)).*(x(2:end)-x(1:end-1)));        
    end
 

