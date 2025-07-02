function varargout = spts(obj,flag1,val1)
 % Computes the Timescales by measuring the "spead" of the system  
 % dE/dt = J*E at E(0) = w(T,Y)
% scale them by F
% (c) 2012 U.P. for VIRTUHCON

% epsilon = 0 <=> all time scales counts
% can be overwritten by parameter
 
% set the default species selection strategy

epsilon = 1e-2;
indexTimes = 1:length(obj.times);
obj.lambda = [];

if isempty(obj.massFractions)&&isempty(obj.initialMassFractions)
     MException('batchReactor:fullTimeScale:emptyMassFractions',...
        ['The mass fraction property is empty',...
        ' you shold solve the reactor first before',...
        ' trying to compute time scales.']).throw;
    
    
elseif isempty(obj.massFractions)
%     fprintf('SPTS: Use initial mass fractions\n');
    % For static computations we use intiial mass fractions
    
        
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure);                
    set(obj.gas,'Y',obj.initialMassFractions); 
     
    [DF,F] = obj.sensitivitY([obj.initialTemperature;obj.initialMassFractions]);         
    nf = norm(F,2);
    temptau = DF*F; 

    if  norm(temptau)<1e-16||nf<1e-16
        Gamma = NaN;
    else         
         
        tau = 1./abs(temptau);
        Gamma =  nf/norm(temptau);
    end
    
else
    tau = NaN(obj.nSpecies,length(indexTimes));
    Gamma = NaN(1,length(indexTimes));
    for k = 1:length(indexTimes)
        %   

        y = obj.massFractions(:,indexTimes(k));
        [DF,F] = obj.sensitivitY([obj.temp(indexTimes(k));y]);         
        nf = norm(F,2);
        temptau = DF*F; 

        if  norm(temptau)<1e-6||nf<1e-6
            Gamma(k) = NaN;
        else         
            indx = find(abs(F)/max(abs(F))>epsilon);
            tau(1:length(indx),k) = 1./abs(temptau(indx));
            Gamma(k) =  nf/norm(temptau);
        end
    end   
end 

 

% loop over the time index

obj.lambda =  Gamma;
 

% output
switch nargout
    case 0
        % doin' nothing  
    case 1
        varargout{1} = obj;
    case 2
        varargout{1} = tau;
        varargout{2} = Gamma;
    
    otherwise
        % 
        ME = MException('batchReactor:timeScaleJacobi:WrongNumberOutput',...
            'The number of output arguments must be lesser than three.');
        throw(ME);
end
end
 
 



