function  xing(obj) 
% g.xing crosses the genes of the currentGeneration of genetic information
switch obj.aimsToReduce
    case 'species'
        obj.currentGeneration = obj.genotype.Xing(obj.currentGeneration);
    case 'reactions'
        obj.currentGeneration = obj.genotype.XingR(obj.currentGeneration);
    otherwise
        error('firstGeneration:WrongAim',...
            'The reduction strategie can only be ''species'' or ''reactions''.')
end
end