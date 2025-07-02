function varargout = pvts(obj,species)
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
    indx = obj.mainSpecies(min(1,size(obj.moleFractions,1)),'mole');
elseif ~isempty(obj.massFractions)
    indx = obj.mainSpecies(min(1,size(obj.massFractions,1)),'mass');     
else
     throw(ME);
end

 
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
                indx = obj.getSpeciesIndexByName(species);
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
    
    
     
    domega = zeros(length(indx),length(indexTimes)-1);
    yk = obj.massFractions(:,indexTimes(1:end-1));
    yk1 = yk;
    yk1(indx,:) = obj.massFractions(indx,indexTimes(2:end));    
    for k = 1:length(indexTimes)-1
        % reconstruct the source term   
        
        [~,omega1] = obj.sensitivitY([obj.temp(indexTimes(k));yk1(:,k)]);
        [~,omega]  = obj.sensitivitY([obj.temp(indexTimes(k));yk(:,k)]);
        domega(:,k) = omega1(indx,:)-omega(indx,:);       
    end
     
    dy =(obj.massFractions(indx,indexTimes(2:end))...
        -obj.massFractions(indx,indexTimes(1:end-1))).^2;
    tau =  sqrt(sum(dy)./sum(domega.^2,1)); 
    obj.lambda = [ tau(1) tau];    
end

% output

switch nargout
    case 0
    case 1
        varargout{1} =  tau;    
    otherwise
        % 
end
end
 
 



