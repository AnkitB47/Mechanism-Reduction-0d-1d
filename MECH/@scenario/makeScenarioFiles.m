function makeScenarioFiles(obj,varargin)
% self.makeScenarioFiles()
% 


switch nargin
    case 1 % only self (obj)
        error('makeScenarioFiles:NotEnoughInputArguments',...
            'Not enough input arguments: At least, a chemistry object is needed to read simulation time etc.')
    case 2 %
        if isa(varargin{1},'chemistry')
            chem = varargin{1};
        else
            error('makeScenarioFiles:secondArgument','The first argument must be a chemistry object.')
        end
    otherwise % 
        if isa(varargin{1},'chemistry')
            chem = varargin{1};
        else
            error('makeScenarioFiles:secondArgument','The second argument must be a chemistry object.')
        end
        for k = 3:nargin
            % argument checking
            % hmmm, what arguments...thinking later about...
        end
end
% able to check the reactor
if isempty(chem.reactor)
    % we assume that we make a DRG style file
    fprintf('Make scenario file for DRG\n')
    makeScenarios(obj,chem);
else
    % testing the reactors
    switch chem.reactor
        case {'commonReactor' 'constVolumeReactor' 'constPressReactor'  'counterFlowFlame'}
            error('makeScenarioFiles:notApplicable',...
                'This function  makes DRG style files for SCENGEN, but the reactor is Cantera based.')
        
        otherwise
            % fallback: DRG file style, pure matlab WITHOUT calling scengen
            makeScenarios(obj,chem);
    end
end
end