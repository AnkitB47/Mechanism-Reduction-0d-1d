function val = dist5(obj1,obj2)
%dist4(obj1,obj2)
% jet another dist function

if abs(obj1.ignt-obj2.ignt)>0.25*obj1.ignt
    val = inf;
else
    time = grefine(obj1.times);
    temp1 = interp1(obj1.times,obj1.temp,time);
    temp2 = interp1(obj2.times,obj2.temp,time);

    val = trapez(time,(temp1-temp2).^2)/trapez(time,temp1.^2)...
            +(obj1.ignt-obj2.ignt)^2/obj1.ignt...
                +(obj2.nLivingSpecies/obj1.nLivingSpecies)^2; 
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