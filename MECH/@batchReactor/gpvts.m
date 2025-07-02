function varargout = gpvts(obj,species)
% Computes the progress variable time scale defined in Ihme/Pitsch
% obj = obj.timeScaleProgressVariable(species)
% species must be a cell array containig the names of the species defining
% the progress variable, e.g. the species of the main reaction.
% Example:
%  obj= obj.timesScaleProgressVariable({'CH4' 'O2'});
%  obj= obj.timesScaleProgressVariable() uses the  main species, i.e. the
%  species with the largest value of mole fraction.
% (c) 2012 U.P. for VIRTUHCON

% epsilon = 0 <=> all time scales counts
% can be overwritten by parameter
 
% set the default species selection strategy

ME = MException('batchReactor:timeScaleProgressVariable:emptyMassFractions',...
        ['The mass fraction property is empty',...
        ' you shold solve the reactor first before',...
        ' trying to compute time scales.']); 

 
if ~isempty(obj.moleFractions)
    indx = obj.mainSpecies(4,'mole');
elseif ~isempty(obj.massFractions)
    indx = obj.mainSpecies(4,'mass');
else
     throw(ME);
end
progress = zeros(2,length(obj.times));
indexTimes = 1:length(obj.times);
obj.lambda = [];

if isempty(obj.massFractions)
    throw(ME);
else   
    switch nargin
        case 1
            % without additional parameters
        case 2
           try
                indx =     obj.getSpeciesIndexByName(species)        ;
           catch ME
                throwAsCaller(ME);
            end
        otherwise
            % Hmmmm????
            error('batchReactor:eigenvalueTimeScales:TooManyArguments',...
                'The number of arguments must be 0, 2, or 4.');
    end
    
%     Gamma = NaN(1,length(indexTimes));
    % loop over the time index
    
    c = zeros(1,length(indexTimes));
    omega_c =  zeros(1,length(indexTimes));
    for k = 1:length(indexTimes)
        %        
        y = obj.massFractions(:,indexTimes(k));
        [~,F] = obj.sensitivitY([obj.temp(indexTimes(k));y]);   
        
        c(k) = sum(y(indx));
        omega_c(k) = sum(F(indx));
    end   
end
 
% tau = abs((c(2:end)-c(1))./(omega_c(2:end)-omega_c(1:end-1))); 
tau = abs((c(1:end)-c(1))./omega_c); 
obj.lambda = tau;
% output
switch nargout
    case 0
        figure      
        loglog(obj.times,tau,'m.');
        hold on

    case 1
        varargout{1} =  obj;
    case 2        
        varargout{1} = progress;
        varargout{2} = tau;
    otherwise
        % 
end
end
 
 



