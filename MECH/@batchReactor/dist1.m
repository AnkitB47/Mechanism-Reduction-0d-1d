function val = dist1(obj1,obj2)
%val = dist1(obj1,obj2)
% Simple distance function of two chemistry objects.
% It measures the relative differences in temperature  and in the main species given
% by the inital set in the L2 norm, ie.
% sqrt(sum(int(X_i(obj1)-X_i(obj2))^2))/int(X_i(obj1))
% Note that the norm is scaled by int(obj1). 

% speciesList = obj1.getInitalSpecies;

 
nsp = 6;
[m] = max(obj1.moleFractions,[],2); 
[~,indx] = sort(m,1,'descend');
speciesList = obj1.speciesNames(indx(1:nsp));
fprintf('batchReactor:dist1: Compute error for species ')

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
        val = val+abs(obj1.nSpecies-obj2.nSpecies)/(obj1.nSpecies+obj2.nSpecies);
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
        
        val = sqrt(trapez(time,(temp1-temp2).^2))/sqrt(trapez(time,temp1.^2)) ;
        for k = 1:nsp
            if int(k)<inf
                val = val+int(k);
            end
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