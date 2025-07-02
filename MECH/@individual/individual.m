classdef individual
%individual: class for individuals of a mechanismn polulation
% properties:
%        chromosome: logical
%         reactions: logical
%             fixed: logical
%         mechanism: string
%             value: double
%             names: cell 
%             index: double
% speciesNotInclude: double
% 
% Methods for class individual:
% 
% deleteSpecies   individual      mutationFixedN  xing            
% fix             mutation        removeSpecies   

    
    properties              
        fixed = []; 
        value = inf;     % field for the value of the evaluation
    end
    
    properties (SetAccess = private)
        chromosome = [];   
        reactions = []; % multiplier for reactions, switch on/off      
        index = []; % multiplier for species, switch on/off    
        speciesNotInclude = []; % ~indx to identify removed
                                % species in a fast way
        
    end
    
    % Change properties only by the contructor
    properties (SetAccess = immutable)         
        mechanism = []; % path to mechamism file
        names = [];     % species names       
    end
    
    methods
        function obj = individual(mechfile)
            switch nargin
                case 0
                    % do nothing, empty indvidual object,
                    % only to create (empty) lists etc.
                case 1
                    % constructor method for individual class
                    try
                    obj.names = chemistry.getSpeciesFromMech(mechfile);
                    catch ME
                        throwAsCaller(ME)
                    end
                    isCTI = strcmp(mechfile(strfind(mechfile,'.')+1:end),'cti');
                    if isCTI
                        b = chemistry.indexCtiFile(mechfile);
                        obj.reactions = logical(1:b.nreactions);
                    end
                    obj.chromosome = logical(1:length(obj.names));
                    obj.fixed = ~logical(1:length(obj.names));
                    obj.mechanism = mechfile;                    
                    obj.index = find(obj.chromosome); 
                    obj.speciesNotInclude = find(~obj.chromosome);
                otherwise
            end
            obj.value = inf;
        end
            
        function obj = setMechanism(~,mech)
            switch nargin               
                case 2
                    try
                        obj = individual(mech);
                    catch ME
                        throwAsCaller(ME)
                    end
                otherwise
                    error('individual:setBaseMech','Use setBaseMech setBaseMech(''mechfile'')');
            end            
        end
        function obj = removeSpecies(obj,listOfSpecies)
            % deletes only if it is not a fixed species
            % two possible ways to define the to-remove species:
            for ll = 1:length(obj)
                if isa(listOfSpecies,'numeric')
                    for k = 1:length(listOfSpecies) 
                        if ~obj(ll).fixed(listOfSpecies(k))
                            obj(ll).chromosome(listOfSpecies(k)) = false;
                        end
                    end
                elseif isa(listOfSpecies,'cell')
                    for l = 1:length(obj(ll).names)
                        for k = 1:length(listOfSpecies)                         
                            if ~obj(ll).fixed(l)&&strcmp(obj(ll).names{l},listOfSpecies{k})
                                obj(ll).chromosome(l) = false;
                            end
                        end
                    end
                elseif isa(listOfSpecies,'char')
                    for l = 1:length(obj.names),
                        if ~obj(ll).fixed(l)&&strcmp(obj(ll).names{l},listOfSpecies)
                            obj(ll).chromosome(l) = false;
                        end
                    end
                else
                    %
                end
                obj.index = find(obj.chromosome);
                obj.speciesNotInclude = find(~obj.chromosome);
            end
        end
        
        function obj = deleteSpecies(obj,listOfSpecies)
            % delete action: remove the species also if it is a fixed one
            % two possible ways to define the to-delete species:
            for ll = 1:length(obj)
                if isa(listOfSpecies,'numeric')
                    for k = 1:length(listOfSpecies)  
                        obj(ll).chromosome(listOfSpecies(k)) = false; 
                        obj(ll).fixed(listOfSpecies(k)) = false;                   
                    end
                elseif isa(listOfSpecies,'cell')
                    for l = 1:length(obj(ll).names),
                        for k = 1:length(listOfSpecies)                         
                            if strcmp(obj(ll).names{l},listOfSpecies{k})
                                obj(ll).chromosome(l) = false;
                                obj(ll).fixed(l) = false;
                            end
                        end
                    end
                elseif isa(listOfSpecies,'char')
                    for l = 1:length(obj(ll).names),
                        if strcmp(obj(ll).names{l},listOfSpecies)
                            obj(ll).chromosome(l) = false;
                            obj(ll).fixed(l) = false;
                        end
                    end
                else
                    %
                end
                obj(ll).index = find(obj(ll).chromosome);
                obj(ll).speciesNotInclude = find(~obj(ll).chromosome);
            end
        end
        
        function obj = fix(obj,listOfSpecies)
            % sets species to non-selectable for removing
            if isa(listOfSpecies,'numeric')
                for l = 1:length(obj)
                    for k = 1:length(listOfSpecies)  
                        obj(l).fixed(listOfSpecies(k)) = true;                   
                    end
                end
            elseif isa(listOfSpecies,'cell')
                for ll = 1:length(obj)
                    for l = 1:length(obj(ll).names),                    
                        for k = 1:length(listOfSpecies)                         
                            if strcmp(obj(ll).names{l},listOfSpecies{k})                             
                                obj(ll).fixed(l) = true;
                            end
                        end
                    end
                end
            elseif isa(listOfSpecies,'char')
                for ll = 1:length(obj)
                    for l = 1:length(obj(ll).names),                    
                        if strcmp(obj.names{l},listOfSpecies)                        
                            obj.fixed(l) = true;
                        end
                    end
                end                
            else
                %
            end
        end
        
        function obj = unfix(obj,listOfSpecies)
            for ll = 1:length(obj)
                if isa(listOfSpecies,'numeric')
                    for k = 1:length(listOfSpecies)  
                        obj(ll).fixed(listOfSpecies(k)) = false;                   
                    end
                elseif isa(listOfSpecies,'cell')
                    for l = 1:length(obj(ll).names),
                        for k = 1:length(listOfSpecies)                         
                            if strcmp(obj(ll).names{l},listOfSpecies{k})                             
                                obj(ll).fixed(l) = false;
                            end
                        end
                    end
                elseif isa(listOfSpecies,'char')
                    for l = 1:length(obj(ll).names),
                        if strcmp(obj(ll).names{l},listOfSpecies)                        
                            obj(ll).fixed(l) = false;
                        end
                    end
                else
                    %
                end
            end
        end        
        
        function obj = mutation(obj,prob)
           if nargin <2
               prob = 0.3;
           end
           for ll = 1:length(obj)
               obj(ll).chromosome(~obj(ll).fixed) = mut(obj(ll).chromosome(~obj(ll).fixed),prob);
               obj(ll).index = find(obj(ll).chromosome); 
               obj(ll).speciesNotInclude = find(~obj(ll).chromosome);           
               obj(ll).value = inf;
           end
        end
        
        function obj = mutationN(obj,N)
            for ll = 1:length(obj)
                %obj = obj.mutationN(N)
                % mutates only N entries in obj.chromosome and only if
                % it is not fixed
                indx = 1:length(obj(ll).chromosome);
                for k = 1:N                
                    i = randi(length(indx));
                    if ~obj(ll).fixed(indx(i))
                        obj(ll).chromosome(indx(i)) = ~obj(ll).chromosome(indx(i));
                        indx(i) = [];
                    end
                end
                obj(ll).index = find(obj(ll).chromosome);
                obj(ll).speciesNotInclude = find(~obj(ll).chromosome); 
                obj(ll).value = inf;  
            end
        end
        
        function obj = mutationPlus(obj,N)
            %obj = mutationPlus(obj,N) 
            % Adds randomly N species to the mechanism
            % It loops over all members of obj ...
            for ll = 1:length(obj)
                if isempty(N)
                    N = fix(length(obj(ll).index)/10+1); % approx 10% will be added
                end
                %obj = obj.mutationN(N)
                % mutates only N entries in obj.chromosome and only if
                % it is not fixed
                indx = 1:length(obj(ll).chromosome);
                for k = 1:N                
                    i = randi(length(indx));
                    if ~obj(ll).fixed(indx(i))
                        obj(ll).chromosome(indx(i)) = true;
                        indx(i) = [];
                    end
                end
                obj(ll).index = find(obj(ll).chromosome);
                obj(ll).speciesNotInclude = find(~obj(ll).chromosome); 
                obj(ll).value = inf;  
            end
        end
        
        function obj = mutationMinus(obj,N)
            %obj = mutationMinus(obj,N)
            % Removes randomly N species from the mechanism            
            % It loops over all members  of obj ...
            for ll = 1:length(obj)
                if isempty(N)
                    N = fix(length(obj(ll).index)/10+1); % approx 10% will be added
                end
                %obj = obj.mutationN(N)
                % mutates only N entries in obj.chromosome and only if
                % it is not fixed
                indx = 1:length(obj(ll).chromosome);
                for k = 1:N                
                    i = randi(length(indx));
                    if ~obj(ll).fixed(indx(i))
                        obj(ll).chromosome(indx(i)) = false;
                        indx(i) = [];
                    end
                end
                obj(ll).index = find(obj(ll).chromosome);
                obj(ll).speciesNotInclude = find(~obj(ll).chromosome); 
                obj(ll).value = inf;  
            end
        end
        
        
        function obj = mutationFixedN(obj,N,rounds)
            indx = find(obj.fixed); % indices of fixed species
            free = 1:length(obj.chromosome);
            free(indx) = []; %#ok<FNDSB> % clear the fixed from the mutation list
            nOfSpecies = sum(obj.chromosome); % how many species 
            % are in the mechanism
            nFree = length(free); % how many species can be mutated
            if sum(obj.fixed)>N
                error('mutationWithFixedSpeciesNumber:WrongNumberOfSpecies',...
                    'The number of species in a mechnism must be larger than the number of fixed species')
            end
            if nOfSpecies>N
                % too many species in the mechanism, only deleting
                % allowed.             
                while nOfSpecies>N                
                    cur = randi(nFree) ;
                    obj.chromosome(free(cur)) = false;
                    nOfSpecies = sum(obj.chromosome);
                    obj.index = find(obj.chromosome); 
                    obj.speciesNotInclude = find(~obj.chromosome);
                    obj.value = inf;
                end
            else
                for l = 1:rounds
                    cur = randi(nFree); % a candidate to mutate               
                    if ~obj.chromosome(free(cur))&&(N>nOfSpecies)
                        obj.chromosome(free(cur)) = true;
                        nOfSpecies = sum(obj.chromosome);
                    else
                        obj.chromosome(free(cur)) = false;
                        nOfSpecies = sum(obj.chromosome);
                        %                     fprintf('1=>0\n')    ;            
                    end
                    obj.index = find(obj.chromosome);
                    obj.speciesNotInclude = find(~obj.chromosome);
                    obj.value = inf;                   
                end
            end
            
        end

        function [obj1,obj2] = xing(obj1,obj2)
            %crossing of two individuals
            
            if ~strcmp(obj1.mechanism,obj2.mechanism)
                error('xing:BaseMechError','Both individuals must base on the same mechanism.')
            end
            if length(obj1)>1||length(obj2)>1
                error('individual:xing:OnlyScalarInput','The arguments must be single individual objects.')
            end
            pos1 = randi(length(obj1.chromosome));
            pos2 = randi(length(obj1.chromosome)); 
%             pos2 =  length(obj1.chromosome);
            part1.c = obj1.chromosome(min(pos1,pos2):max(pos1,pos2));
            part1.f = obj1.fixed(min(pos1,pos2):max(pos1,pos2));
            obj1.chromosome(min(pos1,pos2):max(pos1,pos2)) = obj2.chromosome(min(pos1,pos2):max(pos1,pos2));
            obj1.fixed(min(pos1,pos2):max(pos1,pos2)) = obj2.fixed(min(pos1,pos2):max(pos1,pos2));
            obj1.index = find(obj1.chromosome); 
            obj1.speciesNotInclude = find(~obj1.chromosome);
            obj2.chromosome(min(pos1,pos2):max(pos1,pos2)) =  part1.c;
            obj2.fixed(min(pos1,pos2):max(pos1,pos2)) =  part1.f;
            obj2.index = find(obj2.chromosome); 
            obj2.speciesNotInclude = find(~obj2.chromosome);            
            obj1.value = inf;            
            obj2.value = inf;
        end
        
        % methods to manipulate the "reaction chromosome"
        function obj = mutationR(obj,prob)
            % mutates the reaction cromosom
            if nargin <2
                prob = 0.3;
            end
            for ll = 1:length(obj)
                obj(ll).reactions = mut(obj(ll).reactions,prob); 
                obj(ll).value = inf;
            end            
        end
        
        function [obj1,obj2] = xingR(obj1,obj2)
            % x-ing of the reactions cromosom
            
            if length(obj1)>1||length(obj2)>1
                error('individual:xing:OnlyScalarInput','The arguments must be single individual objects.')
            end
            pos1 = randi(length(obj1.reactions));
            pos2 = randi(length(obj1.reactions)); 
            excange = obj1.chromosome(min(pos1,pos2):max(pos1,pos2));
            obj1.reactions(min(pos1,pos2):max(pos1,pos2)) = obj2.reactions(min(pos1,pos2):max(pos1,pos2));
            obj2.reactions(min(pos1,pos2):max(pos1,pos2)) =  excange;
        end
        
        
    end
    
end

