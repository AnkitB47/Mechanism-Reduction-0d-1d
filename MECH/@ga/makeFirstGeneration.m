function  makeFirstGeneration(obj,nind)
% g.firstGeneration creates a first generation by mutation of the
% base individual. It fills the currentGeneration properties
% with a vector of indviduals.

if nargin < 2
    nind = 2*length(obj.genotype.chromosome);
elseif isempty(nind)
    nind = 2*length(obj.genotype.chromosome);
end
fprintf(['First generation consists of ',...
    num2str(nind),' individuals\n']);
if ~isa(obj.genotype,'individual')
    error('Xing:InPutMustBeIndividual',...
        'The input must be a cell array. Each entry must be a individual object')
end
% call the individual mutation method
for k = 1:nind
    switch obj.aimsToReduce
        case 'species'
            obj.currentGeneration = ...
                [obj.currentGeneration obj.genotype.mutation(obj.options.prob)];  
        case 'reactions'
            obj.currentGeneration = ...
                [obj.currentGeneration obj.genotype.mutationR(obj.options.prob)];
        otherwise
            error('firstGeneration:WrongAim',...
                'The reduction strategie can only be ''species'' or ''reactions''.')
    end
end
end