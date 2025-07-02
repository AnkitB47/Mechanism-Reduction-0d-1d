function val = dist3(obj1,obj2)
%val = dist1(obj1,obj2)
% Simple distance function of two chemistry objects.
% It measures the relative differences btw. obj1 and obj2  only in the temperature.

 

time = grefine(obj1.times);


temp1 = interp1(obj1.times,obj1.temp,time);
temp2 = interp1(obj2.times,obj2.temp,time);

val = trapez(time,(temp1-temp2).^2)/trapez(time,temp1.^2)...            
                +(obj2.nSpecies/obj1.nSpecies)^2; 


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