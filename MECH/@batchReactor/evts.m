function varargout = evts(obj,flag1,val1,flag2,val2)
% Eigenvalue based Time Scale inspired by Jean Caudal...
% (c) 2012 U.P. for VIRTUHCON
% 
% br.timeScaleCaudal(['eps',val,'step',val])
% Usage: 
% Input parameter
% a) empty
% b) 'eps',val
%  uses val as "crtiteria of relevance" to select time scales 
% c) 'steps',val 
%   evaluates every val-th time step or evaluates at times val if
%   length(val)>1
% Examples
%  br.timeScaleCaudal('eps',0.01);
%  br.timeScaleCaudal('eps',0.01,'steps',2);
%  br.timeScaleCaudal(['eps',1e-3,'steps',[0 0.1 0.5 1])
% Output parameter
% a) empty plots relevant time scales
% b) one gives back the batchReactor-object with obj.lambda = Tau
% c) two [tau,Tau]  where tau are all and Tau are the relevant time scales 
% Examples
%  br.timeScaleCaudal('eps',0.01); 
% plot relevant time scales
%  br = br.timeScaleCaudal('eps',0.01); 
% updates the object
%  [tau,Tau] = br.timeScaleCaudal('eps',0.01);
%  loglog(br.times,tau,'y.')
%  hold on
%  loglog(br.times,Tau,'k.')
% computes time scales and plots all in yellow and relevant in black 
% Note eps  = 0 <=> all time scales are relevant
 

epsilon = 0.01;
indexTimes = 1:length(obj.times);
obj.lambda = [];
if isempty(obj.massFractions)&&isempty(obj.initialMassFractions)
    ME = MException('batchReactor:fullTimeScale:emptyMassFractions',...
        ['The mass fraction property is empty',...
        ' you shold solve the reactor first before',...
        ' trying to compute time scales.']);
    throw(ME);
    
elseif isempty(obj.massFractions)
    fprintf('SPTS: Use initial mass fractions\n');
    % For static computations we use intiial mass fractions
    obj.massFractions = obj.initialMassFractions;
    obj.temp = obj.initialTemperature;
    indexTimes = 1;
end   
    switch nargin
        case 1
            % without additional parameters
        case 3
            switch flag1
                case 'steps'
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
            end
        case 5
            switch flag1
                case 'steps'
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
            end
            switch flag2
                case 'steps'
                    if length(val2)>1
                        indexTimes = getTimesIndex(obj,N );
                    else
                        indexTimes = 1:val2:length(obj.times);
                        fprintf(['Take every ',num2str(val2),'-th timestep, ',...
                            'i.e. evaluate ',num2str(length(indexTimes)),' points\n'])
                    end
                 case {'eps' 'epsilon'}
                    epsilon = val2;
                otherwise
            end            
        otherwise
            % Hmmmm????
            ME = MException('batchReactor:timeScaleCaudal:TooManyArguments',...
                'The number of arguments must be 0, 2, or 4.'); 
            throw(ME);
    end
    tau = zeros(obj.nSpecies,length(indexTimes));   
%     gamma = zeros(obj.nSpecies,length(indexTimes));   
     
    Tau = NaN(obj.nSpecies,length( indexTimes));    
 
    for k = 1:length(indexTimes)
        %        
        [DF,F] = obj.sensitivitY([obj.temp(indexTimes(k));obj.massFractions(:,indexTimes(k))]);    
        [V,eigenval] = eig(DF,'nobalance');          
        l = 1;
        while l <=obj.nSpecies
             if isreal(eigenval(l,l))
                 tau(l,k) = 1/abs(eigenval(l,l));  
                 % leave V as it was and jumpt to the next EV            
                 l = l+1;                
             else
                 % Causals complex rule : 
                 % tau(i) = 1/| Re(lambda)| and 
                 % tau(i+1) = 1/| Im(lambda) |
                 tau(l,k) = 1/abs(real(eigenval(l,l)));                
                 tau(l+1,k) = 1/abs(imag(eigenval(l,l))); 
                 % put the real and imaginary part of the Eigenvector in
                 % the new basis V
                 V(:,l:l+1) = [real(V(:,l)),imag(V(:,l))];                 
                 % jump over the conjugated brother of lambda
                 l = l+2;
             end             
        end
        % norm the eigenvectors becaus 'no balance" option gives back an
        % non normed eigenvector basis.
        V = V*diag(sqrt(1./diag(V'*V)));        
        % V is now real  and normal
        
        % compute the coefficients of  omega aka F in basis V  
        % because of V normal gamma = V\F
        gamma = V\F;

        % scaling, caudal's version...        
        normgamma =  abs(gamma)/max(abs(gamma));        

        % select gamma  wrt epsilon 
        indx = normgamma>=epsilon; 
        Tau(indx,k) = tau(indx,k);  
        %       relative approximation error
        % e(k) = norm(F-V(:,indx)*gamma(indx))/norm(F);
    end        
    obj.lambda = Tau;
 
% output
switch nargout
    case 0
        % do nothing 
    case 1
        % object output
        varargout{1} = obj;
    case 2   
        % tau/Tau output
        varargout{1} = tau;
        varargout{2} = Tau;
    otherwise
        % 
end
end





