function val = dist1(obj1,obj2)
% obj1.dist1(obj2) measures the distance btw. two flame1d objects
% wrt. the mole fractions of the main species and the temperature

% compue the "main" species

nsp = 3;
[m] = max(obj1.massFractions,[],2); 
[~,indx] = sort(m,1,'descend');
speciesList = obj1.speciesNames(indx(1:nsp));
fprintf('flame1d:dist1: Compute error for species ')
for k = 1:nsp-1
    fprintf([speciesList{k},', ']);
end
fprintf([speciesList{k+1},'\n']);

    
space = grefine(obj1.space); 
n = length(speciesList) ;
for k = 1:nsp           
    moleFrac1 = interp1(obj1.space,obj1.getMoleFractionsByName(speciesList{k}),space);
    moleFrac2 = interp1(obj2.space,obj2.getMoleFractionsByName(speciesList{k}),space);
    int(k) = 1/nsp*sqrt(trapez(space,(moleFrac1-moleFrac2).^2))/sqrt(trapez(space,moleFrac1.^2)); 
end
temp1 = interp1(obj1.space,obj1.temp,space);
temp2 = interp1(obj2.space,obj2.temp,space);

 

val = sqrt(trapez(space,(temp1-temp2).^2))/sqrt(trapez(space,temp1.^2)) ;
% % % for k = 1:n
% % %     if int(k)<inf
% % %         val = val+int(k);
% % %     end
% % % end
figure(10)
plot(space,[temp1;temp2]);drawnow
fprintf('Difference in temperature %d\n',val)
if val > 0.1
    val = inf;
else
    val = obj2.nSpecies+val;
end

function grid = grefine(grid)
        % refines the discrete space
        grid = sort([grid 0.5*(grid(2:end)+grid(1:end-1))]);
    end

    function val = trapez(x,y)
        % trapezidual rule
         val = sum(0.5*(y(2:end)+y(1:end-1)).*(x(2:end)-x(1:end-1)));        
    end
end