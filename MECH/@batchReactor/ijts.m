function varargout = ijts(obj,flag1,val1)
% Old time scale definition 
% Computes the Timescales by approximating the EV by the diagonal of DF and
% scale them by F
% (c) 2012 U.P. for VIRTUHCON

% epsilon = 0 <=> all time scales counts
% can be overwritten by parameter
 
% set the default species selection strategy

epsilon = 1e-2;
indexTimes = 1:length(obj.times);
obj.lambda = [];

if isempty(obj.massFractions)
    error('batchReactor:fullTimeScale:emptyMassFractions',...
        ['The mass fraction property is empty',...
        ' you shold solve the reactor first before',...
        ' trying to compute time scales.'])
else   
    switch nargin
        case 1
            % without additional parameters
        case 3
            switch flag1
                case 'step'
                    if length(val1)>1
                        indexTimes = getTimesIndex(obj,val1 );
                    else
                        indexTimes = 1:val1:length(obj.times);
                        fprintf(['Take every ',num2str(val1),'-th timestep, ',...
                            'i.e. evaluate ',num2str(length(indexTimes)),' points\n'])
                    end   
                case {'eps' 'epsilon'}
                    epsilon = val1;
                   
                otherwise
                   error('batchReactor:eigenvalueTimeScales:UnkownProperty',...
                        ['The property ',flag1,' is unknown.' ]);  
            end            
        otherwise
            % Hmmmm????
            error('batchReactor:eigenvalueTimeScales:TooManyArguments',...
                'The number of arguments must be 0, 2, or 4.');
    end
    Gamma = NaN(obj.nSpecies,length(obj.times));
    tau = NaN(obj.nSpecies,length(obj.times));
%     Gamma = NaN(1,length(indexTimes));
    % loop over the time index
    for k = 1:length(indexTimes)
        %        
        y = obj.massFractions(:,indexTimes(k));
        [DF,F] = obj.sensitivitY([obj.temp(indexTimes(k));y]);          
        d = abs(diag(DF));     
        indx = find(abs(F)/max(abs(F))>=epsilon);                
        Gamma(1:length(indx),k) = 1./abs(d(indx)); 
        tau(:,k) = abs(1./d);       
    end   
    obj.lambda = abs(Gamma);
end

% output
switch nargout
    case 0
% do nothing
    case 1
        varargout{1} = obj;
    case 2        
        varargout{1} = tau;
        varargout{2} = obj.lambda;
    otherwise
        % 
end
end
 
 



