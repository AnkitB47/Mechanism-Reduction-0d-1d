function createMech(obj,baseMech,mechFile)
% Method for creating a mechanism file from the
% genetic information of individual class
%  Ususage:
% individual.createMech(baseMech,mechFile)
% baseMech should be the same mechanism as in 
% individual.mechanism, but it may differ in the format cti/inp or the file
% location. mechFile ist the name of the traget. It can be a valid path.
% (C) 2012 Uwe Prüfert for VIRTUHCON

% only to instiate a son of chemistry
reactor = batchReactor(baseMech);
try
    reactor.removeSpecies(baseMech,mechFile,...
                        obj.speciesNotInclude,...
                        obj.names(obj.speciesNotInclude));
catch ME
    throw(ME);
end
end