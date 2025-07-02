function   solve(obj,varargin)
% obj.solve(['reactorType'])
% solve the batch reactor defined in obj.reactor or via the optional parameter.
% If none ist given, commonReactor will be called
% Examples:
%     r1.reactor = 'constVolumeReactor';
%     r1 = solve  
% solves the contant volume bathc reactor.
%     r1 = solve('tibox')
% solves the constant pressure reactor by the tibox programm.  
switch nargin
    case 1
        reactorToCall = obj.reactor; 
        if isempty(reactorToCall)
            reactorToCall = ' ';                  
        end
    case 2
        reactorToCall = varargin{1};
    otherwise
        ME = MException('batchReactor:solve:WrongNumberInput',...
            ['Wrong number of input arguments.'...
            ' Only one additional argument is suppored.']);
        throwAsCaller(ME);
        
end
%
switch reactorToCall
    case 'constVolumeReactor'
        obj.constVolumeReactor;
    case 'constVolumeNetReactor'
        obj.constVolumeNetReactor;
    case 'constPressReactor'
        obj.constPressReactor;
    case 'constPressNetReactor' 
        obj.constPressNetReactor;  
    case 'constTempReactor' 
        obj.constTempReactor;      
    case 'tibox'                    
        obj.tibox;
    otherwise
        obj.commonReactor;
end
% When solve the reactor, the time scales are not longer valid
obj.lambda = [];
end