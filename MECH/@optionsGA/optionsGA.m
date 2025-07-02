classdef optionsGA
    %DATA class for GA options
    %    
    
    properties
        ngen;
        selectionStrategie;
        selectionThreshold;
        nind;  
        prob;
    end
    
    methods
        function obj = optionsGA()
            obj.ngen = 5;
            obj.selectionStrategie = 'survivors';
            obj.selectionThreshold = 6;
            obj.nind = 50;   
            obj.prob = 0.3;
        end
    end
    
end

