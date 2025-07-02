function species = selectSpeciesFromList(elementsLine)
% utility function for removeSpecies
% it reads the species from the elementsLine and gives back a cell array
% containing the species names


l = 0;
inList = false;
spec = [];
for k = 1:length(elementsLine)
    
    
    if inList && ~(strcmp(elementsLine(k),' ')||strcmp(elementsLine(k),'"'))
        spec = [spec,elementsLine(k)];
    end
    if inList && (strcmp(elementsLine(k),' ')||strcmp(elementsLine(k),'"'))
        l = l+1;
        if ~isempty(spec)
            candidate{l}= spec;
            spec = [];
        end
    end 
    if strcmp(elementsLine(k),'"')
       inList = swtch(inList);
    end
end


l = 1;
for k = 1:length(candidate)
    if ~(isempty(candidate{k})||strcmp(candidate{k},'=>')||strcmp(candidate{k},'<=>')||strcmp(candidate{k},'+'))
        species{l}=candidate{k};
        l = l+1;
    end
end




    function b = swtch(b)       
        b = ~b;                
    end
end