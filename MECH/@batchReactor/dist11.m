function val = dist11(obj1,obj2)
%val = dist1(obj1,obj2)
% Distance function of two chemistry objects.
% It measures the relative differences in temperature  and in the main
% species. If the error in the species and temperature is larger than a
% tolerance of 20%, val = Inf, else val = sum (int(Y_i-Y_i^mech)) #species^mech / #species_ref
 

 
% nsp = 6;
% [m] = max(obj1.moleFractions,[],2); 
% [~,indx] = sort(m,1,'descend');
[~,speciesList] = obj1.mainSpecies;
nsp = length(speciesList);
fprintf('\nCompute error for species:\n')

for k = 1:nsp-1
    fprintf([speciesList{k},', ']);
end
fprintf([speciesList{k+1},'\n']);

switch obj1.reactor    
    case {'constVolumeNetReactor' 'constPressNetReactor'}
        % reactors with constant number of grid points
        time =  obj1.times;
        for k = 1:nsp 
            int(k) = 1/nsp*sqrt(trapez(time,(obj1.getMoleFractionsByName(speciesList{k})...
                -obj2.getMoleFractionsByName(speciesList{k})).^2))...
                /sqrt(trapez(time,obj1.getMoleFractionsByName(speciesList{k}).^2)); 
        end
        
        val = sqrt(trapez(time,(obj1.temp-obj2.temp).^2))/sqrt(trapez(time,obj1.temp.^2)) ;
        for k = 1:nsp
            if int(k)<inf
                val = val+int(k);
            end
        end
        fprintf(['Error in species = ',num2str(val),'\n']);
        if val<1
            val = val+1-abs(obj1.nSpecies-obj2.nSpecies)/(obj1.nSpecies+obj2.nSpecies);
        else
            val = inf;
        end         
    otherwise
        % reactors with non-constant number of grid points.
        % use a refined mesh for both objects
        time = grefine(obj1.times);       
        for k = 1:nsp            
            moleFrac1 = interp1(obj1.times,obj1.getMoleFractionsByName(speciesList{k}),time);
            moleFrac2 = interp1(obj2.times,obj2.getMoleFractionsByName(speciesList{k}),time);
            int(k) = 1/nsp*sqrt(trapez(time,(moleFrac1-moleFrac2).^2))/sqrt(trapez(time,moleFrac1.^2)); 
        end
        temp1 = interp1(obj1.times,obj1.temp,time);
        temp2 = interp1(obj2.times,obj2.temp,time);
        val = sum(int);
       
        fprintf(['Error in species = ',num2str(val),'\n']);
        if val <1
            val = val+1-abs(obj1.nSpecies-obj2.nSpecies)/(obj1.nSpecies+obj2.nSpecies)+...
                sqrt(trapez(time,(temp1-temp2).^2))/sqrt(trapez(time,temp1.^2)) ;
        else
            val = inf;
        end
end

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
