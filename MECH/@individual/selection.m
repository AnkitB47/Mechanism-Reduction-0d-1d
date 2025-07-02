function individualsListOut = selection(~,individualsListIn,varargin)
%individualsListOut = selection(individualsListIn,[argument,value]) 
% selection takes care that no identical indviduals  appear in the output
% list 
switch length(varargin)
    case 0
        numberOfSurvivors = min(4,length(individualsListIn));
        cases = 2;
    case 2
        if ~isnumeric(varargin{2})
            error('selection:ArgumentMustBeNumeric','The second argument must be a scalar')
        else
            if strcmp(varargin{1},'threshold') 
                cases = 1;
            elseif strcmp(varargin{1},'survivors')
                cases = 2;
                numberOfSurvivors = min(length(individualsListIn),varargin{2});
                if length(individualsListIn)<varargin{2}
                    error('selection:ArgumentFitting',...
                        'The number of survivors is larger than the number of members in the population.')
                end
            else               
                error('selection:ArgumentIsNotValid',...
                    'The argument is not defined: It must be eigther threshold or survivors')
            end
        end        
    otherwise
        %
        error('selection:WrongNumberOfArguments','There are more than two or only one argument.')
end
        
n = length(individualsListIn);
switch cases
    case 1 % all individuals survive iff their value of fitness is <= threshold 
        m = 0;
        for k = 1:n           
            if individualsListIn(k).value<=varargin{2}
                unique = true;
                mm = m;
                for l = 1:mm                     
                    if individualsListIn(k).chromosome==individualsListOut(l).chromosome
                        unique = false; 
                        break
                    end
                end
                if unique
                   m = m+1;
                   individualsListOut(m) = individualsListIn(k);
                end
            end
        end                
    case 2 % we select the numberOfSurvivors-th fittest
        values = zeros(1,n);
        for k = 1:n     % copy the .value in a  vector
            values(k) = individualsListIn(k).value;
        end
        [~,indx] = sort(values);
        m = 0;
        for k = 1:n            
            unique = true;
            mm = m;
            for l = 1:mm
                if individualsListIn(indx(k)).chromosome==individualsListOut(l).chromosome
                    unique = false; 
                    break
                end
            end
            if unique
                m = m+1;
                individualsListOut(m) = individualsListIn(indx(k));
            end
            if m>=numberOfSurvivors 
                break
            end
        end   
    otherwise
        error('individual:selction:StrangeError',...
            'Ooops, strange things happen here ... you NEVER-EVER  should see this message!')
end
end
