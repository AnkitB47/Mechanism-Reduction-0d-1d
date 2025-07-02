function val = dist13(obj1,obj2,varargin)
%val = dist13(obj1,obj2)
% Computes the distance between two batchReactor objects. 
% implements the error(y-y) < 0.1 condition


% define the weightf for the measures:
% alpha [species temperature lambda #species]
if nargin==2
    alpha = [1,0,0,1];
else
    alpha = varargin{1};
end



[~,speciesList] = obj1.mainSpecies(4);
nsp = length(speciesList);
 
if isempty(speciesList)&& alpha(1)>0
    error('batchReactor:dist13:EmptyFractions','The Species list is empty.');
elseif alpha(1)>0
    fprintf('batchReactor:dist13: Compute error for species:\n ');
    for k = 1:nsp-1
        fprintf([speciesList{k},', ']);
    end
    fprintf([speciesList{k+1},'\n']);
end





switch obj1.reactor    
    case {'constTempReactor', 'constPressReactor','constVolumeReactor'}
        % reactors with non-constant number of grid points.
        % use a refined mesh for both objects
        time = unique([obj1.times  obj2.times]);
        
        [~,M,~] = assema1d(time,[1 1 1]',0,1,0);
        M = diag(sum(M));
        int = zeros(1,nsp+2); % species + temp + lambda
        if  alpha(1)>0 %only for alpha ~= 0
            for k = 1:nsp        
                moleFrac1 = interp1(obj1.times,obj1.getMoleFractionsByName(speciesList{k}),time);
                moleFrac2 = interp1(obj2.times,obj2.getMoleFractionsByName(speciesList{k}),time);     
                if sqrt(moleFrac1*M*moleFrac1')<1e-12
                    warning('batchReactor:dist13:IntegralAlmostZero',...
                        'The Integral is almost zero, value of relative error may be inacurate.')
                end
                err = sqrt((moleFrac1-moleFrac2)*M*(moleFrac1-moleFrac2)')/...
                        sqrt(moleFrac1*M*moleFrac1');
                    fprintf(['dist13: Species ',speciesList{k},' has relative error of ',num2str(err*100),' %%\n']);
                if err<0.2
                    int(k) = 0; 
                else
                    int(k) = Inf;           
                end
            end
        end
        % compute the timescales if neccessary 
        if isempty(obj1.lambda)&&alpha(3)~=0
           fprintf('dist13:Compute Timescales for input one\n');
           obj1 = obj1.eigenvalueTimeScales('eps',1e-3);
        end      
        if isempty(obj2.lambda)&&alpha(3)~=0
            fprintf('dist13:Compute Timescales for input two\n');
            obj2 = obj2.eigenvalueTimeScales('eps',1e-3);
        end
        
        if alpha(2)>0
            temp1 = interp1(obj1.times,obj1.temp,time);
            temp2 = interp1(obj2.times,obj2.temp,time);        
            if  ((temp1-temp2)*M*(temp1-temp2)')/...
                    (temp1*M*temp1')<0.2
                int(nsp+1) = 0;
            else
                int(nsp+1) = Inf;
            end;
        end   
         
        
        if alpha(3)>0
            timescale1 = interp1(obj1.times,obj1.lambda,time);
            timescale2 = interp1(obj2.times,obj2.lambda,time);
            int(nsp+2) =  alpha(3) * sqrt((timescale1 -  timescale2)*M*(timescale1 -  timescale2)');
        end
        
        
        
        val = sum(int);
    otherwise
        error('batchReactor:WrongReactorType','Wrong reactor type for this measure.')
end
val = val + alpha(4)*obj2.nSpecies/obj1.nSpecies;
 

     
   
end