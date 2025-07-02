function solveForReactions(obj)
% Method for reducing mechanism by GA on the reactions. The reactions are
% switched on/off by using the setMultiplier method of kinetics class.

% coppy the reference object
%current = obj.reference;
% re-define the reactions
obj.aimsToReduce = 'reactions';
if isempty(obj.currentGeneration)
    obj.makeFirstGeneration(250);
end
obj.evaluate; 
end