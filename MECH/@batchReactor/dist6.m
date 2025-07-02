function val = dist6(obj1,obj2)
%dist4(obj1,obj2)
% jet another dist function

if abs(obj1.ignt-obj2.ignt)/obj1.ignt>0.33 || abs(obj1.temp(end)-obj2.temp(end))/obj1.temp(end) > 0.1
    val = inf;
else    
    val = (obj2.nSpecies/obj1.nSpecies)^2+obj2.nSpecies/obj2.nLivingSpecies; 
end
end