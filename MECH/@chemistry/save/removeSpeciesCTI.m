function removeSpeciesCTI(mechFileIn,mechFileOut,deadSpeciesIndexList,deadSpeciesNamesList)
% removeDeadSpecies(mechFileIn,mechFileOut,deadSpeciesIndexList)
%   Removes the species (and Reactions) defined in deadSpeciesIndexList and
%   deadSpeciesNamesList from mechFileIn, the resulting mechanism will be
%   written in mechFileOut. 
%   The input arguments deadSpeciesIndexList,deadSpeciesNamesList must be
%   consistent, i.e. it holds  
%          species(deadSpeciesIndexList) = deadSpeciesNamesList.



% check for stupid inputs/output file declaration
if strcmp(mechFileIn,mechFileOut)
    error('removeSpeciesCTI:SameInputOutputFile','The input and the output files are the same, are you sure? However, it is not allowed. Bye!')
end
% try to open the files

fprintf('removeSpeciesCTI: Indexing input file\n')
try     
    indexCti = chemistry.indexCtiFile(mechFileIn);
catch ME
    throw(ME)     
end
fprintf(['removeSpeciesCTI: Mechanism File contains ', num2str(length(indexCti.species)),' species ',...
    num2str(length(indexCti.reactions)+length(indexCti.three_bodies)+length(indexCti.falloffs)),' reactions.\n'])


% check the input list
if ~isnumeric(deadSpeciesIndexList)
    error('removeSpeciesCTI:ListMustBeNumeric','Dead Species Index List must contain the INDICES of dead species')
else
    % list is ok
end
fprintf(['removeSpeciesCTI: There are ',num2str(length(deadSpeciesIndexList)),' species in dead species list.\n'])
fprintf('removeSpeciesCTI: I found \n\t')
listTheSpecies(deadSpeciesNamesList)
fprintf('removeSpeciesCTI: I will now try to delete all these species and related reacion.\n\t It will take only a few of your earth seconds.\n')


%----------------------------------------------------
% some logical field to mark the "areas of action"

% specify the lines we have to work or not to work
% initialize all true
nothingToDo = logical(1:indexCti.nlines);
handleReactions = logical(1:indexCti.nlines);   handleReactions = ~handleReactions;
handleIdealGas = logical(1:indexCti.nlines);    handleIdealGas = ~handleIdealGas;
notToCopy = logical(1:indexCti.nlines); notToCopy = ~notToCopy;




% but we have to modify all ideal_gas blocks, mark the "to removed species"
nothingToDo(indexCti.ideal_gas(1):indexCti.species(1)-2) = false;
% we have an action for the to remove species..
for k = 1:length(deadSpeciesIndexList)
    if deadSpeciesIndexList(k)==length(indexCti.species)
        % the last species is choosen, then we mark all lines until the
        % reactions block starts 
        nothingToDo(indexCti.species(deadSpeciesIndexList(k)):...
            min([indexCti.falloffs(1),indexCti.three_bodies(1),indexCti.reactions(1)])-8) = false;
    else    % the regular case
        nothingToDo(indexCti.species(deadSpeciesIndexList(k)):indexCti.species(deadSpeciesIndexList(k)+1)-1) = false;
    end
end

% and for all reactions
nothingToDo(indexCti.reactions) = false; % mark reactions be a canditate   
nothingToDo(indexCti.reactions+1) = false; % mark the options line, it can be empty                                           
nothingToDo(indexCti.three_bodies) = false; % mark reaction to  be a canditate 
nothingToDo(indexCti.three_bodies+1) = false; % mark the efficiencies line                                                                                              
nothingToDo(indexCti.falloffs) = false;       % mark falloffs to be a canditate
nothingToDo(indexCti.falloffs+1) = false;     % there are four line after the falloff_reaction keyword
nothingToDo(indexCti.falloffs+2) = false;     % we have to test
nothingToDo(indexCti.falloffs+3) = false;
nothingToDo(indexCti.falloffs+4) = false; 
handleIdealGas(indexCti.ideal_gas(1):indexCti.species(1)-2) = true;

% the special actions:
% now we remove all species defined in  the parameter "deadSpeciesIndexList"
% thew simplest action: to handle the "to remove species" 
for k = 1:length(deadSpeciesIndexList)
    if deadSpeciesIndexList(k)==length(indexCti.species)
        % the lst species is choosen, then we mark all lines until the
        % reactions block starts
        notToCopy(indexCti.species(deadSpeciesIndexList(k)):...
            min([indexCti.falloffs(1),indexCti.three_bodies(1),indexCti.reactions(1)])-8) = true;
    else    % the regular case
        notToCopy(indexCti.species(deadSpeciesIndexList(k)):indexCti.species(deadSpeciesIndexList(k)+1)-1) = true;
    end
end
handleReactions(indexCti.reactions) = true;
handleReactions(indexCti.three_bodies) = true;
handleReactions(indexCti.falloffs) = true;
inSpecies = false;
deadReactionRemove = false;
threeBodyRemove = false;
falloffRemove  = false;

% open the files
[fidMech,mesg] = fopen(mechFileIn,'r');
if fidMech < 0 
    fprintf(mesg)
end 
[fidNewMech,mesg] = fopen(mechFileOut,'w');
if fidNewMech < 0 
    fprintf(mesg)
end

% for counting the actions
nToRemoveReactions = 0;
% a counter to make the ideal_gas mames unique
nGasBlocks = 0;

for k = 1:indexCti.nlines
    if nothingToDo(k)   
        %-----------------------------------
        % comments living species etc
        % just copy...
        %-----------------------------------
        cline = fgetl(fidMech); 
        fprintf(fidNewMech,[cline,'\n']);
        %
    elseif handleIdealGas(k)
        %-----------------------------------
        % the ideal gas block contains also a species block that should be
        % modified. We copy all lines exept this species block. However, its
        % length is not fixed so we must scan it
        % The species block starts with the keyword "species" and ends with the
        % keyword "reactions". The keywords start at pos 7 in every line
        cline = fgetl(fidMech); 
        if length(cline)>=13 && strcmp(cline(7:13),'species')
            inSpecies = true; % switch species block "on"
        elseif length(cline)>=15&& strcmp(cline(7:15),'reactions')
            inSpecies = false; % switch species block "off"
        end
        if inSpecies
            % the "species action"
            % analyse the line and remove the dead species
            nline = chemistry.removeDeadSpeciesFromIdealGasBlock(cline,deadSpeciesNamesList);
            fprintf(fidNewMech,[nline,'\n']);            
        elseif length(cline)>=13 && strcmp(cline(7:14),'kinetics')
            % can be omitted                %
        elseif length(cline)>=13 && strcmp(cline(7:15),'transport')           % 
            fprintf(fidNewMech,[cline,'\n']);   
        elseif length(cline)>=13 && strcmp(cline(1:5),'ideal')
            nGasBlocks = nGasBlocks + 1;
            fprintf(fidNewMech,['ideal_gas(name = "gas_',num2str(nGasBlocks),'",\n']);              
        else
            fprintf(fidNewMech,[cline,'\n']);
        end            
    elseif notToCopy(k)    
        %-----------------------------------
        % The  definition lines of dead species:
        % removing by "comment out"
        %-----------------------------------
%         cline = fgetl(fidMech);
        fgetl(fidMech); 
%         fprintf(fidNewMech,['#  ',cline,'\n']);  
    elseif handleReactions(k)
        %-----------------------------------
        % The definition lines  of dead recations: 
        % removing  by "comment out" 
        %-----------------------------------
        cline = fgetl(fidMech);
        
        % ckecking if the reaction contains a dead species
        if chemistry.isDeadReaction(cline,deadSpeciesNamesList)
%             fprintf(fidNewMech,['#  ',cline,'\n']);
            deadReactionRemove = true;
            nToRemoveReactions = nToRemoveReactions +1;
%            
            if strcmp(cline(1:5),'three') 
                % mark the line after a dead three body as "to be removed"
                threeBodyRemove = true;
            elseif strcmp(cline(1:4),'fall') 
                falloffRemove = true;
            end
        else
            % if the reaction contains only living species, we copy it into
            % the reduced mech file
            deadReactionRemove = false;
            falloffRemove = false;
            threeBodyRemove = false;
            fprintf(fidNewMech,[cline,'\n']);
        end
    else
        cline = fgetl(fidMech);
        if length(cline)>16&&strcmp(cline(10:16),'options')&&deadReactionRemove
%             fprintf(fidNewMech,['#',cline,'\n']);  
            deadReactionRemove = false;
        elseif threeBodyRemove
            % remove the efficnecies line of a three body reaction
            threeBodyRemove = false;
%             fprintf(fidNewMech,['#',cline,'\n']);           
        elseif falloffRemove
            % remove the falloff lines 2-5
            % analyse the first two characters in the lines
            if length(cline)>=11
                key = cline(10:11);
            else
                key = ' ';
            end
            switch key
                case {'kf' 'fa'}
%                     fprintf(fidNewMech,['#',cline,'\n']);                    
                case 'ef'
%                     fprintf(fidNewMech,['#',cline,'\n']);   
                    falloffRemove = false;
                otherwise
            end
        elseif length(cline)>=21,             
            if strcmp(cline(10:21),'efficiencies') 
                % remove the dead species from the efficiencies list
                newEfficiencies = chemistry.removeDeadFromEfficiencies(cline,deadSpeciesNamesList);
                fprintf(fidNewMech,[newEfficiencies,'\n']);
            else
                fprintf(fidNewMech,[cline,'\n']);
            end
        else
            fprintf(fidNewMech,[cline,'\n']);
        end
        % up to now we write a empty line instead to make the "real action"
        % here we define here the action    
    end
end
fprintf('removeDeadSpecies: Try to write the files and close them.\n')
% close all files

st = fclose(fidNewMech);
if st<0
    error('removespeciesCTI:fcloseFails','fclose fals to close the mechfile')
end
st = fclose(fidMech);
if st<0
    error('removespeciesCTI:fcloseFails','fclose fals to close the mechfile')
end
st = fclose('all');
if st<0
    error('removespeciesCTI:fcloseFails','fclose fals to close the mechfile')
end
 




fprintf(['removeDeadSpecies: Afterall, I remove ',num2str(nToRemoveReactions),...
    ' of ',num2str(length(indexCti.reactions)+length(indexCti.three_bodies)+length(indexCti.falloffs)),' reactions'])
fprintf(['...or round ',num2str(round(100*nToRemoveReactions/(length(indexCti.reactions)+length(indexCti.three_bodies)+length(indexCti.falloffs)))),' percent. \n'])
fprintf('removeDeadSpecies: OK. Thats all, see you!\n\n')  
end


% a tiny helper function: it lists a cell array full of strings
function listTheSpecies(list)
    m = ceil(length(list)/4);
    r = length(list);
    n = 1;
    fprintf('\n')
    for k = 1:m,
        for l = 1:4,
            if n <= r
                fprintf(['\t',list{n}]);
                n = n+1;
            end
        end
        fprintf('\n')
    end
    fprintf('\n')
end
