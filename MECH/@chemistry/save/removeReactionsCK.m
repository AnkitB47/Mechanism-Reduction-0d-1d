function removeReactionsCK(mechFileIn,mechFileOut,varargin)
% removeReactionsCK(mechFileIn,mechFileOut[,species_list])
% remove all reactions that contains the species in species_list.
% species_list must be a CELL array, e.g. {'H2' 'CO2'} or a species object.
% If no species_list is given, removeReactions removes ALL reactions from 
% the mechanism.
% (c) U. Pruefert (VIRTUHCON) 2011.

% check for stupid inputs/output file declaration
if strcmp(mechFileIn,mechFileOut)
    error('removeDeadSpecies:SameInputOutputFile','The input and the output files are the same, are you sure? However, it is not allowed. Bye!')
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

% some comment add to the mechanism file
fprintf(fidNewMech,['! This is a reduced mechanism based on ',mechFileIn,' \n']);
fprintf(fidNewMech,'! The removing process is performed by removeReactionsCK\n');
fprintf(fidNewMech,'! \t\t (c) U. Pruefert 2011.\n!\n');

if isempty(varargin)
    all = true;
else
    if isa(varargin{1},'cell')
        deadSpeciesNamesList = varargin{1}; 
    elseif isa(varargin{1},'species')
        deadSpeciesNamesList = varargin{1}.names;
    else
        error('removeReactionsCK:wrongArgument','The optinal argument must be a CELL array containing the names of species or a SPECIES object')
    end
    all = false;
    fprintf(fidNewMech,'! All reactions associated with the following species will be deleted:\n');
    listTheSpecies(fidNewMech,deadSpeciesNamesList)
end
inReactions = false;
toRemove = true;
nOfReactions = 0;nOfRemovedReactions = 0;
while 1 % loop over all lines
    cline = fgetl(fidMech); 
    if ~ischar(cline),         
        break,
    end 
    if (length(cline)>=4) && (strcmp(cline(1:4),'REAC')||strcmp(cline(1:4),'reac'))
        % reactions block
        inReactions = true;
        fprintf(fidNewMech,[cline,'\n']);
    elseif (length(cline)>=3) && (strcmp(cline(1:3),'END')||strcmp(cline(1:3),'end')) 
        % exit all blocks         
        inReactions = false;
        
        fprintf(fidNewMech,[cline,'\n']);
    else
        % now, we have to work...
        % Definitely, we are in a non keyword line or a (regular) comment
        % line
        if inReactions
            % first we look for a line that defines a reaction.
            % if a reaction contains a to-removed-species, all lines can be
            % removed untill we found the next reaction
            if isReactionLine(cline)
                if ~all
                    toRemove = false;
                end
                nOfReactions = nOfReactions + 1;
            end
            if ~toRemove
                if containsDeadSpecies(cline,deadSpeciesNamesList)
                    toRemove = true;
                    nOfRemovedReactions = nOfRemovedReactions + 1;
                end
            end
            if ~toRemove                
                fprintf(fidNewMech,[cline,'\n']);
            end
        else
            fprintf(fidNewMech,[cline,'\n']);            
        end
    end   
end

% we close all files
err = fclose(fidNewMech);
if err < 0
    fprintf('An error during closing the file is occured. Please check the file names etc.\n')
end
err = fclose(fidMech);
if err < 0
    fprintf('An error during closing the file is occured. Please check the file names etc.\n')
end


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
        % local local functions for identifying species by name
        function b = isLeft(cline,l)
            % checks if a species string is a whole species name or only a part
            % of a species name: left of the name...
            b = (l==1) || ((l>1) && (strcmp(cline(l-1),'+') ||  strcmp(cline(l-1),'>')||...
                strcmp(cline(l-1),'=') || (double(cline(l-1))>48 && double(cline(l-1))< 58) )...
                 &&...
                ((l>=2)||(strcmp(cline(l-2),'+') ||  strcmp(cline(l-2),'>') || strcmp(cline(l-2),'=')))); 
        end
        % of a species name: right of the name...
        function b = isRight(cline,l)
            b = ( strcmp(cline(l+lngthSpecName),' ') ||  strcmp(cline(l+lngthSpecName),'<') ||...
                strcmp(cline(l+lngthSpecName),'+')||  strcmp(cline(l+lngthSpecName),'=') ||...
                strcmp(cline(l+lngthSpecName:l+lngthSpecName+1),'(+')) ; 
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
end