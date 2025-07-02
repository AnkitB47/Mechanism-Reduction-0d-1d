function varargout = timeScaleCaudalMod(obj,flag1,val1)
% Eigenvalue based Time Scale inspired by Jean Caudal...
% In contrast to the original method the mod function uses a different
% relevance criteria.
% (c) 2012 U.P. for VIRTUHCON
% 
% br.timeScaleCaudalMod(['N',val])
% Usage: 
% Input parameter
% a) empty
% b) 'N',val
%  uses val as "crtiteria of relevance" to select N time scales 
%  After computing the order of reducing the error ||w-sum(a_i*v_i)||
% Examples
%  br.timeScaleCaudal('N',3);
% 
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
obj.lambda = [];
N = 3;
indexTimes = 1:length(obj.times);

if isempty(obj.massFractions)
    ME = MException('batchReactor:fullTimeScale:emptyMassFractions',...
        ['The mass fraction property is empty',...
        ' you shold solve the reactor first before',...
        ' trying to compute time scales.']);
    throw(ME);
else   
    switch nargin
        case 1
            % without additional parameters
        case 3
            switch flag1                
                case {'N' 'n'}
                    N = val1;
                otherwise
                    ME = MException('batchReactor:timeScaleCaudalMod:UnkownArgument',...
                        ['The argument ''',flag1,''' is unknown.']); 
                    throw(ME);
            end
        otherwise
            % Hmmmm????
            ME = MException('batchReactor:timeScaleCaudalMod:TooManyArguments',...
                'The number of arguments must be 0 or 2.'); 
            throw(ME);
    end
    tau = zeros(obj.nSpecies,length(indexTimes));   
%     gamma = zeros(obj.nSpecies,length(indexTimes));   
     
    Tau = NaN(1,length(indexTimes));    
   
    for k = 1:length(indexTimes)
        %        
        [DF,F] = obj.sensitivitY([obj.temp(indexTimes(k));obj.massFractions(:,indexTimes(k))]);    
        [V,eigenval] = eig(DF,'nobalance');    
 
        l = 1;
        while l <=obj.nSpecies
             if isreal(eigenval(l,l))
                 tau(l,k) = 1/abs(eigenval(l,l));                
                 l = l+1;
             else
                 % complex rule : 
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
         
        % select gamma  wrt epsilon 
        indx = greedy_search(F,gamma,V);   
        Tau(1:N,k)= tau(indx(1:N),k);           
    end 
    
    obj.lambda = Tau;
   
    
end
% output
switch nargout
    case 0
        % do nothing but plot the timescales
        loglog(obj.times,Tau,'.')
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

function indx = greedy_search(omega,a,v)
% min (sum (i)) 
% st.omega - sum a_i*v_i
% function that gives a "optimal" index vector for the
% approximation of omega by a linear combination of the
% sum a_i*v_i

% compute a_i *v_i
AV = v*diag(a);
na = length(a);
y = 0;
indx = [];
coeffNr = 1:na;
abs_minimum = Inf;
for l = 1:na
    aa = zeros(1,length(coeffNr));     
    for k = 1:length(coeffNr)          
        aa(k) = norm(omega-(y+AV(:,k)))/norm(omega);     
    end
    % compute the largest descent
    [minimum,kk] = min(aa);
    if minimum>abs_minimum
       break
    else
        abs_minimum = minimum;
    end
    y = y+AV(:,kk);
    indx = [indx coeffNr(kk)];    
    AV(:,kk) = []; 
    coeffNr(kk) = [];
end
end



