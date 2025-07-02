function  mutate(obj,prob)
% mutates the genetic information
if nargin<2
    prob = 0.3;
end
switch obj.aimsToReduce
    case 'species'
        for k = 1:length(obj.currentGeneration)
            obj.currentGeneration(k) = ...
                obj.currentGeneration(k).mutation(prob);  
        end
    case 'reactions'
        % not jet implemented, but comming soon
        for k = 1:length(obj.currentGeneration)
            obj.currentGeneration(k) =...
                obj.currentGeneration(k).mutationR(prob);  
        end
    otherwise
        error('firstGeneration:WrongAim',...
            'The reduction strategie can only be ''species'' or ''reactions''.')
end
end