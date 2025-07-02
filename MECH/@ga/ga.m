 classdef ga < handle
%   ga  class for genetic algorithm to minimize 
%   min f(x)
%   st. g(x) = 0
%   where g(x) solves a chemical reactor, flame etc. based on chemical
%   source term reduction.

%   properties
    properties(Access=public)
        genotype@individual = individual(); 
        objective@function_handle;
        reference@chemistry;
        aimsToReduce@char;
        options@optionsGA = optionsGA(); 
    end
    
    properties(SetAccess=private)
        currentGeneration; % do NOT initialize here anything!               
    end
    
    properties (Constant,Access = private)
        % define some errors
         wrInpFile = MException('ga:WrongInputFileFormat',...
                                    'The input must be a cti, inp or xml File.');
         wrNumbInp = MException('ga:WrongNumberInputs',...
                        'The number of input arguments must be one.')
         noString = MException('ga:NoString',...
                            'The argument must be a string');       
         noFile = MException('ga:NoFileFound',...
                            'Cannot open the mechanism file.');
    end
    
%   
%   methods 


    methods(Access = public)
        % hight level (user) methods: constructor, solve...
        function obj = ga(varargin)
            % constructor method for ga
            switch length(varargin)
                case 0
                    obj.genotype = individual();
                    obj.objective = '';  
                    obj.aimsToReduce = 'species';               
                case {1 2 3}
                    for k = 1:length(varargin)
                        switch class(varargin{k})
                            case 'individual'
                                obj.genotype = varargin{k};
                            case 'function_handle'
                                obj.objective = varargin{k};
                            case {'batchReactor' 'flame1d' 'reactorNet'}
                            obj.reference = varargin{k};
                            otherwise
                                ME = MException('ga:WrongClassOfArguments',...
                                    ['Arguments must be of class ''individual'','...
                                    ' ''function_handle'' or ''batchReactor''.']);
                                throw(ME);                        
                        end
                        obj.aimsToReduce = 'species';  
                    end                    
                otherwise
                    throw(obj.wrNumbInp)                    
            end                
        end        
        solve(obj,options);
        solveForReactions(obj);
        guessFixSpecies(obj,nsp);
    end  
 
    methods(Access = private)
        % low level methods, mainly used by solve    
        makeFirstGeneration(obj,nind); 
        mutate(obj,prob);
        mutateOne(obj)  
        xing(obj);
        select(obj,varargin);
        evaluate(obj);
    end
end
