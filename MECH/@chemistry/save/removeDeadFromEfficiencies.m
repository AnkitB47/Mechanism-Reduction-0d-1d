function newEfficiencies = removeDeadFromEfficiencies(oldEfficiencies,deads)
%
% scans efficiencies list and remove species in dead
l = 1;m=1;remove = false;
startName = false;
for k = 26:length(oldEfficiencies)
    if isLetter(oldEfficiencies(k)) && ~startName
        startSpeciesName(l) = k;
        startName = true;
        l = l+1;
    elseif strcmp(oldEfficiencies(k),':')
        endSpeciesName(m) = k-1;
        startName = false;
        m = m+1;
    end
end
endSpecies = startSpeciesName(2:end)-3;
endSpecies(length(startSpeciesName))=length(oldEfficiencies)-3;

if ~(length(startSpeciesName)==length(endSpeciesName))
    MException('chemistry:WrongLengthArguments',...
        'The format of the input file may be corrupted.').throw;
end
newEfficiencies = oldEfficiencies(1:25) ;
for k = 1:length(startSpeciesName),
    for l = 1:length(deads)  
        if strcmp(oldEfficiencies(startSpeciesName(k):endSpeciesName(k)),deads{l}),
            remove = true;           
            break
        else
            remove = false;
        end
    end
    if ~remove
        newEfficiencies = [newEfficiencies,' ',oldEfficiencies(startSpeciesName(k):endSpecies(k)),' ' ];
    end
end
newEfficiencies = [newEfficiencies,'")']  ;
function b = isLetter(charSpec)
    % is it a  letter character (to define a species name, the first must
    % be a letter, 20H is NOT a well formed species name
    b = (double(charSpec)<91 && double(charSpec)>64) ||...% capitals A...Z
                (double(charSpec)<123 && double(charSpec)>96);  % lowercasses a..z
end    
end

