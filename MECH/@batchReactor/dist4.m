function val = dist4(obj1,obj2)
%dist4(obj1,obj2)
% jet another dist function

if abs(obj1.ignt-obj2.ignt)>0.25*obj1.ignt
    val = inf;
else
    val = (obj1.ignt-obj2.ignt)^2+(obj2.nLivingSpecies/obj1.nLivingSpecies); 
end



end

