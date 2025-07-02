function M = sensitivityMatrix(obj,t)
indx = obj.getTimesIndex(t);
M = obj.sensitivitY([obj.temp(indx);obj.massFractions(:,indx)]);
end