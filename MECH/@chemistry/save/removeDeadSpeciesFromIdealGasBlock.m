function nline = removeDeadSpeciesFromIdealGasBlock(reactionLine,deads)
%nline = removeDeadSpeciesFromIdealGasBlock(reactionLine,deads)
%   We remove the species from the dead list and write back a new species
%   definition line

% extract the species from the line after the keyword
species = chemistry.scanSpeciesList(reactionLine(20:end));

newSpecies=[];
 
 
for l = 1:length(species)
    for m = 1:length(deads)
        if strcmp(species{l},deads{m})
           species{l} = '';         
        else
            %fprintf(['remove species ',species{l},' from list\n'])
        end
    end
end
l = 1;
for k = 1:length(species)
    if ~isempty(species{k})
        newSpecies{l}=species{k};
        l = l+1;
    end
end
newSpeciesString = [];
for k = 1:length(newSpecies)
    newSpeciesString = [newSpeciesString,' ',newSpecies{k},' '];
end
% we add the new line to the first 19 entries in the old line, i.e we write
% " species" keyword and preserve the structure of the file...
% we also have to set the block closing '""",'
if strcmp(reactionLine(end-3:end),'""",')
    nline = [reactionLine(1:18),newSpeciesString,'""",'];
elseif strcmp(reactionLine(7:13),'species')
    
    nline = [reactionLine(1:19),newSpeciesString];
else
    
    nline = [reactionLine(1:19),newSpeciesString];
end
end

