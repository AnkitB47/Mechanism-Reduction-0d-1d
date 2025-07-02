function val = distPartikelTemp(obj1,obj2)
 

time = grefine(obj1.times);
time = grefine(time);

           
    
temp1 = interp1(obj1.times,obj1.temp,time);
temp2 = interp1(obj2.times,obj2.temp,time);
 intT1 = trapez(time,(temp1-obj1.initialTemperature).^2);
intErr = trapez(time,(temp1-temp2).^2);
val =  sqrt(intErr/intT1);

fprintf(['\nNorm of temperature = ',num2str(sqrt(intT1)),'\n']);
fprintf(['Intergral error (T-T_mech)^2  = ',num2str(intErr),'\n']);
fprintf(['Relative error in temperature = ',num2str(val),'\n']);
fprintf(['Number of species in mechanism = ',num2str(obj2.nSpecies),'\n\n'])

if val>0.01
    val = Inf;
else
    val = obj2.nSpecies+val;
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
 

