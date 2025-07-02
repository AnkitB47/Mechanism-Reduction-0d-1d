function  select(obj,varargin)
% selects the best ones from the current generation of individuals
switch length(varargin)
    case 0
        switch obj.aimsToReduce
            case 'species'
                obj.currentGeneration = ...
                    obj.genotype.selection(obj.currentGeneration);
            case 'reactions'
                % to implement 
            otherwise
        end
    case 2
        switch obj.aimsToReduce
            case 'species'
                obj.currentGeneration = ...
                    obj.genotype.selection(obj.currentGeneration,...
                    varargin{1},varargin{2});
            case 'reactions'
                % to implement 
            otherwise
        end
    otherwise
        ME =MException('ga:select:WrongNumberOfArguments',...
            'The number of arguments must be 0 or 2.');
        throw(ME);
end
end