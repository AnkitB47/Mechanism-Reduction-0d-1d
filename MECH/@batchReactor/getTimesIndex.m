function index= getTimesIndex(obj,times )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isempty(obj.times)
    error('batchReactor:getTimesIndex:emptyTimesField',...
        'The times property is empty. You should first call solve')
end

if max(times)>max(obj.times)||min(times)<min(obj.times)
    error('batchReactor:getTimesIndex:outOfRange',...
        'One or more values in input times are aut of the times in object.times');
end

index = zeros(size(times));
for k = 1:length(times)
    [~,b] = min(abs(obj.times -times(k)));
    index(k) = b;
end

