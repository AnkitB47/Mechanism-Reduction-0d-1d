function obj = setTimes(obj,times)
% Sets the time vector, empty argument clears the 
% times property
% obj = obj.setTimes(times)
% obj = setTimes()
if nargin==2
    obj.times = times;
else
    obj.times = [];
end
end