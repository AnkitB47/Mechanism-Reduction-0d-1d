classdef species
    %species   class: a container for cell arrays containig species names.
    % The main idea is to have  operators like \  etc.
    %    
    
    properties
        names = [];
        mech = [];
    end
    
    methods
        function obj = species(listORstring)
            % two constructor call possible
            if isa(listORstring,'cell')
                % constructor
                obj.names = listORstring;
            elseif isa(listORstring,'char')
                % input is a mechanism file
                [obj.names,obj.mech] = getSpeciesFromMech(listORstring);
            end
        end
        
        function obj = mldivide(left,right)
            % A  \ B operator                   
            for k = 1:length(left.names)
                for l = 1:length(right.names)
                   if strcmp(left.names{k},right.names{l})
                        left.names{k} ='';
                        
                   end
                end
            end
            l = 1;
            for k = 1:length(left.names)
                if ~isempty(left.names{k})
                    spName{l}=left.names{k};
                    l = l+1;
                end
            end
            obj = species(spName);
        end
        
        function n = size(obj)
            n(1) = 1;
            n(2) = length(obj.names);
        end
    end
    
end

