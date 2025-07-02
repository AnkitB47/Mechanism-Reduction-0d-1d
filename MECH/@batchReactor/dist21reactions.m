function val = dist21reactions(obj1,obj2)
 

time = grefine(obj1.times);
time = grefine(time);

           
    
temp1 = interp1(obj1.times,obj1.temp,time);
temp2 = interp1(obj2.times,obj2.temp,time);
val =  sqrt(trapez(time,(temp1-temp2).^2)/trapez(time,temp1.^2));

fprintf(['Relative error in temperature = ',num2str(val),'\n']);

if val>0.3
    val = Inf;
else
    val = obj2.nActiveReactions;
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
 

