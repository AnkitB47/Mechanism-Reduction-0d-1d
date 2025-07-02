function removeSpeciesCK(mechFileIn,mechFileOut,varargin)
%removeSpeciesCK(mechFileIn,mechFileOut,deadSpeciesNamesList)
% This is the CHEMKIN Version of removeSpecies.
%   Removes the species (and Reactions) defined in 
%   deadSpeciesNamesList from mechFileIn, the resulting mechanism will be
%   written in mechFileOut. 
%   The input arguments deadSpeciesNamesList must be
%   either a CELL array containing the names of the species, e.g. {'OH' 'CH3'}
%   or a species object  
%   (c) U. Pruefert (VIRTUHCON) 2011.
 


% check for stupid inputs/output file declaration
if strcmp(mechFileIn,mechFileOut)
    error('removeDeadSpecies:SameInputOutputFile','The input and the output files are the same, are you sure? However, it is not allowed. Bye!')
end

if isa(varargin{1},'cell')
    deadSpeciesNamesList = varargin{1}; 
elseif isa(varargin{1},'species')
    deadSpeciesNamesList = varargin{1}.names;
else
    error('removeSpeciesCK:wrongArgument','The optinal argument must be a CELL array containing the names of species or a SPECIES object')
end


% try to open the files
[fidMech,mesg] = fopen(mechFileIn,'r');
if fidMech < 0 
    fprintf(mesg)
end 
[fidNewMech,mesg] = fopen(mechFileOut,'w');
if fidNewMech < 0 
    fprintf(mesg)
end

% initalize the line counter
nLines = 0;
nOfReactions = 0;
nOfRemovedReactions = 0;

% some important counters
inElements = false;
inSpecies = false;
inThermo = false;
inReactions = false;
linesJetToRemove = 0;

% write some comments in the mechanism file
fprintf(fidNewMech,['! This is a reduced mechanism based on ',mechFileIn,' \n']);
fprintf(fidNewMech,'! The following species are deleted:\n');
listTheSpecies(fidNewMech,deadSpeciesNamesList)
fprintf(fidNewMech,'! Also, all reactions that contain species listed above are removed.\n');
fprintf(fidNewMech,'! The removing process is performed by removeSpeciesCK\n');
fprintf(fidNewMech,'! \t\t (c) U. Pruefert 2011.\n!\n');



while 1 % loop over all lines
    cline = fgetl(fidMech); 
    % because we use a fixed condition in our while loop, we need a loop
    % break statement: if the next "line"contains NO character, it can only
    % be a EOF
    if ~ischar(cline),         
        break,
    end  
    nLines = nLines + 1;
    if (length(cline)>=1) && strcmp(cline(1),'!')
        % comment line, write it in the output or not???
        fprintf(fidNewMech,[cline,'\n']);
    elseif (length(cline)>=4) && (strcmp(cline(1:4),'ELEM')||strcmp(cline(1:4),'elem')) 
        inElements = true;
        fprintf(fidNewMech,[cline,'\n']);
    elseif (length(cline)>=4) && (strcmp(cline(1:4),'SPEC')||strcmp(cline(1:4),'spec')) 
        % species block: Remove all species from the to-remove-list
        inSpecies = true;
        fprintf(fidNewMech,[cline,'\n']);
    elseif (length(cline)>=4) && (strcmp(cline(1:4),'THER')||strcmp(cline(1:4),'ther')) 
        inThermo = true;
        fprintf(fidNewMech,[cline,'\n']);
    elseif (length(cline)>=4) && (strcmp(cline(1:4),'REAC')||strcmp(cline(1:4),'reac'))
        % reactions block
        inReactions = true;
        fprintf(fidNewMech,[cline,'\n']);
    elseif (length(cline)>=3) && (strcmp(cline(1:3),'END')||strcmp(cline(1:3),'end')) 
        % exit al block
        inSpecies = false;
        inThermo = false;
        inReactions = false;
        inElements = false;
        fprintf(fidNewMech,[cline,'\n']);
    else
        % now, we have to work...
        % Definitely, we are in a non keyword line or a (regular) comment
        % line
        if inElements 
            % only to copy the element block 
            fprintf(fidNewMech,[cline,'\n']);            
        elseif inSpecies
            % remove all species from the line and fill it with blanks
            fprintf(fidNewMech,[removeSpeciesFromLine(cline,deadSpeciesNamesList),'\n']);   
        elseif inThermo
            % a thermo block has exactly 4 lines. If we identify a to-rmove
            %-species in a block, the next three lines are also be removed
            if thermoBlockToRemove(cline,deadSpeciesNamesList)
                linesJetToRemove = 3;
            elseif linesJetToRemove > 0               
                linesJetToRemove = linesJetToRemove - 1;
            else
                fprintf(fidNewMech,[cline,'\n']);
            end
        elseif inReactions
            % first we look for a line that defines a reaction.
            % if a reaction contains a to-removed-species, all lines can be
            % removed untill we found the next reaction
            if isReactionLine(cline)
                toRemove = false;
                nOfReactions = nOfReactions + 1;
            end
            if ~toRemove
                if containsDeadSpecies(cline,deadSpeciesNamesList)
                    toRemove = true;
                    nOfRemovedReactions = nOfRemovedReactions + 1;
                end
            end
            if ~toRemove
                % maybe we must remove species from the efficiencies line
                if isEfficienciesLine(cline) && efficienciesContainsDeadSpecies(cline,deadSpeciesNamesList)
                    cline = removeSpeciesFromEfficiencies(cline,deadSpeciesNamesList);                   
                end
                fprintf(fidNewMech,[cline,'\n']);
            end
        end
    end
    
      
end

% we close all files
fclose('all');

 

% some statitics
fprintf('\tremoveSpeciesCK: I removed the following ') 
fprintf([num2str(length(deadSpeciesNamesList)),' species:\n']);
listTheSpecies(1,deadSpeciesNamesList)
fprintf(['\n and\n ',num2str(nOfRemovedReactions),' of ',num2str(nOfReactions),' reactions.\n']);


% And now: some internal functions
% we replace the to-removed-species in the species block by blanks
    function cline = removeSpeciesFromLine(cline,deadSpeciesNamesList)
        for k = 1:length(deadSpeciesNamesList)
            lngthSpecName = length(deadSpeciesNamesList{k});
            for l = 1:length(cline)-lngthSpecName + 1;                
                foundSpec = strcmp(cline(l:l+lngthSpecName-1),deadSpeciesNamesList{k});
                if (((l+lngthSpecName)>length(cline) ||  strcmp(cline(l+lngthSpecName),' ')) && foundSpec)...
                        && ((l==1 ||  strcmp(cline(l-1),' ')) && foundSpec)
                    cline(l:l+lngthSpecName-1) = blank(lngthSpecName);
                end                 
            end
        end        
    end

% we identify the species in the thermo definition line
    function b = thermoBlockToRemove(cline,deadSpeciesNamesList)
        b = false;
        for k = 1:length(deadSpeciesNamesList)
            lngthSpecName = length(deadSpeciesNamesList{k});
            b = strcmp(cline(1:lngthSpecName),deadSpeciesNamesList{k})&&strcmp(cline(lngthSpecName+1),' ');
            if b               
                return
            end
            
        end         
    end

% we check if a line is a reaction definition line
    function b = isReactionLine(cline)
        b = false;
        for l = 1:length(cline)-1;
            % if we found a comment mark, we can abort the procedure
            if strcmp(cline(l),'!')
                return
            end
            % a reaction line contains "=" (within "=>", "<=", "<=>", or "=" )
            b = strcmp(cline(l),'=');
            if b 
                return                          
            end
        end
    end

% we check, if a line contains a species
    function b = containsDeadSpecies(cline,deadSpeciesNamesList)
        b = false;
        for k = 1:length(deadSpeciesNamesList)
            lngthSpecName = length(deadSpeciesNamesList{k});
            for l = 1:length(cline)-lngthSpecName-1;
                foundSpec = strcmp(cline(l:l+lngthSpecName-1),deadSpeciesNamesList{k});
                b = isLeft(cline,l) && isRight(cline,l) && foundSpec ; 
                if b 
                    return
                end               
            end
        end 
        % local functions for identifying species by name
        function b = isLeft(cline,l)
            % checks if a species string is a whole species name or only a part
            % of a species name: left of the name...
            b = (l==1) || ((l>1) && (strcmp(cline(l-1),'+') ||  strcmp(cline(l-1),'>')||...
                strcmp(cline(l-1),'=') || (double(cline(l-1))>48 && double(cline(l-1))< 58) )...
                 &&...
                ((l>=2)||(strcmp(cline(l-2),'+') ||  strcmp(cline(l-2),'>') || strcmp(cline(l-2),'=')))); 
        end
        % of a species name: right of the name can be a blank, "<", "+" "=" "(+", or "/"...
        function b = isRight(cline,l)
            b = ( strcmp(cline(l+lngthSpecName),' ') ||  strcmp(cline(l+lngthSpecName),'<') ||...
                strcmp(cline(l+lngthSpecName),'+')||  strcmp(cline(l+lngthSpecName),'=') ||...
                strcmp(cline(l+lngthSpecName:l+lngthSpecName+1),'(+') ||...
                strcmp(cline(l+lngthSpecName),'/')) ; 
        end
    end

% helper function that lists the species
    function listTheSpecies(fid,list)
        m = ceil(length(list)/4);
        r = length(list);
        n = 1;
        fprintf(fid,'!\n');
        for k = 1:m
            fprintf(fid,'!');
            for l = 1:4,
                if n <= r
                    fprintf(fid,['\t',list{n}]);
                    n = n+1;
                end
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'! \n');
    end

% looking for efficencies  lines
    function b = isEfficienciesLine(cline)
        % an efficiencies line uses a slash as delimiter 
        b = ~isempty(findstr(cline,'/'));        
    end

% remove a Species from the efficiencies line
    function cline = removeSpeciesFromEfficiencies(cline,deadSpeciesNamesList)
        for k = 1:length(deadSpeciesNamesList)
            startSpecies = min(findstr(cline,[deadSpeciesNamesList{k},'/']));
            if ~isempty(startSpecies)              
                endspecies = min(findstr(cline(startSpecies:end),'/'));
                endspecies = startSpecies+endspecies;
                endspecies = endspecies+min(findstr(cline(endspecies+1:end),'/'));
                cline(startSpecies:endspecies) = blanks(endspecies-startSpecies+1);
            end       
        end
    end
% fill with blanks
    function str = blank(n)
            str = [];
            for kk = 1:n
                str = [str,' '];
            end
    end
    function b =  efficienciesContainsDeadSpecies(cline,deadSpeciesNamesList)
        b = false;
        for k = 1:length(deadSpeciesNamesList)
            lngthSpecName = length(deadSpeciesNamesList{k});
            for l = 1:length(cline)-lngthSpecName-1;
                foundSpec = strcmp(cline(l:l+lngthSpecName-1),deadSpeciesNamesList{k});
                b = isLeft(cline,l) && isRight(cline,l) && foundSpec ; 
                if b 
                    return
                end               
            end
        end 
        % local functions for identifying species by name
        function b = isLeft(cline,l)
            % checks if a species string is a whole species name or only a part
            % of a species name: left of the name...
            b = (l==1) || (strcmp(cline(l-1),' ')) ;
        end
        % of a species name: right of the name can be a blank, "<", "+" "=" "(+", or "/"...
        function b = isRight(cline,l)
            b = ( strcmp(cline(l+lngthSpecName),' ') ||  strcmp(cline(l+lngthSpecName),'<') ||...
                strcmp(cline(l+lngthSpecName),'+')||  strcmp(cline(l+lngthSpecName),'=') ||...
                strcmp(cline(l+lngthSpecName:l+lngthSpecName+1),'(+') ||...
                strcmp(cline(l+lngthSpecName),'/')) ; 
        end
    end
end

