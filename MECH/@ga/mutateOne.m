function  mutateOne(obj)   
% mutates only one bit of the genetic information.
switch obj.aimsToReduce
    case 'species'
        for k = 1:length(obj.currentGeneration)
            obj.currentGeneration(k) = obj.currentGeneration(k).mutationN(1);  
        end
    case 'reactions'
        % not jet implemented, but comming soon
    otherwise
        error('firstGeneration:WrongAim',...
            'The reduction strategie can only be ''species'' or ''reactions''.')
end
end