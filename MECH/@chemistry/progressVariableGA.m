function [val,alpha] = progressVariableGA(obj)
% progressVariableGA - computes a progress variable wrt. the
% massFractions. It uses a Genetic Algorithm with 12 generations and
% 100 individuals.
% Note: This is an ALPHA version of the code.

    if isempty(obj.massFractions)
        MException('chemistry:progressVariableGA:EMPTYMASSFRACTIONS',...
            'Empty mass fractions').throw;
    end
    
    y = obj.massFractions(obj.getLivingSpeciesIndex,:);  
    % 
    fprintf('Create initial population for GA\n');
    [alpha,values] = firstGen(y);  

    fprintf('NoGeneration\t Value of f\n');
    fprintf(['\t',num2str(0),'\t',num2str(values(1)),'\n']);
    fprintf('Start evolution: Survival of the  fittest\n');
    for k = 1:20
        [alpha,values] = darwinism(alpha,values,y);        
        fprintf(['\t',num2str(k),'\t',num2str(values(1)),'\n']);
        if readyOrNot(values,'stopp')
            break
        end
    end

    val = z(alpha(1,:),y);
    alpha = alpha(1,:);
end

% some local stuff

function [alpha,values] = darwinism(alpha,values,y)
% darwinism - runs a generation in the cirlce of life of PVs
% IN  : alpha, y
% OUT : alpha and values of objective
    nSelectedInd = 10;
    ns = size(y,1);
    % dummi line: only to remember...
    % alpha(1:10,nSelectedInd) = alpha(1:10,nSelectedInd);
    % crossing 
    for k = 1:6
        while true % to prevent to have only one point.
            indx1 = unique(randi([1,ns],1,2));
            if length(indx1)>1
                break  
            end
        end
        % three sections crossing
        alpha(k*nSelectedInd+1:(k+1)*nSelectedInd,1:indx1(1)) = ...
            alpha(randperm(nSelectedInd),1:indx1(1));
        alpha(k*nSelectedInd+1:(k+1)*nSelectedInd,indx1(1)+1:indx1(2)) = ...
            alpha(randperm(nSelectedInd),indx1(1)+1:indx1(2));
        alpha(k*nSelectedInd+1:(k+1)*nSelectedInd,indx1(2)+1:end) = ...
            alpha(randperm(nSelectedInd),indx1(2)+1:end);
    end
    for k = 7:9
        % two segment crossing
        indx1 =  randi([1,ns],1,1);      
        alpha(k*nSelectedInd+1:(k+1)*nSelectedInd,1:indx1) = ...
            alpha(randperm(nSelectedInd),1:indx1);        
        alpha(k*nSelectedInd+1:(k+1)*nSelectedInd,indx1+1:end) = ...
            alpha(randperm(nSelectedInd),indx1+1:end);
    end
    for k = 1:250
        % random mutation of the whole bunch
        % Note that here we can produce non-advantaged mutation of the
        % champions, hence the values can be increase in the next
        % generation.
        column = randi([1,ns]);
        row = randi([1,100]);
        val = -1+randi([0,2],1,1);
        alpha(row,column) = val;
    end
    for k = 1:100
        values(k) = objective(alpha(k,:),y); 
    end    
    [ values,indx] = sort(values);
    alpha = alpha(indx,:);  
end

function [alpha,values] = firstGen(y)
% firstGen - generates the inital generation for GA
% IN  : data y
% OUT : alpha, values of objective(alpha,y) 
%       in increasing order.
    nIndFirst = 1000;
    ns = size(y,1);
    alpha = zeros(nIndFirst,ns);
    val = zeros(1,nIndFirst);
    for k = 1:nIndFirst     
        alpha(k,:) = -1+randi([0,2],1,ns);
        val(k) = objective(alpha(k,:),y); 
    end    
    [ values,indx] = sort(val);
    alpha = alpha(indx,:);     
end

function val = z(alpha,y)
% z - normalized progress variable
% IN  : coefficients alpha and data y
% OUT : values z(x)
    val = alpha * y;
    val = (val-min(val));
    val = 1/max(val)*val;
end

function val = dz(alpha,y)
% dz - difference z(k+1)-z/k)
% IN  : coefficients alpha and data y
% OUT : z(k+1)-z(k) for all x(k). Nx-1 values!
   z1 = z(alpha,y);    
   val = z1(2:end)-z1(1:end-1);    
end


function val = DzY(alpha,y)
% DzY - derivative dY/dz
% IN  : coefficients alpha and data y
% OUT : matrix ns x nx-1 dy/dz
    DY = y(:,2:end)-y(:,1:end-1);
    val = DY*diag(1./dz(alpha,y));
end

function val = f(y)
% f - function that handles the value DzY
% IN  : data y
% OUT : value of f. sum(sum(y.^2)).
%     val = sum(sum(y.^2));
    val = sum(max(abs(y)),2);
end

function val = valid(alpha,y)
% valid - function that handles 
%         the monotonicity constraint
% IN  : coefficients alpha and data y
% OUT : zero if it is feasible, Inf if not
    if all(dz(alpha,y)>0)
        val = 0;
    else
        val = Inf;
    end
end

function val = objective(alpha,y)
% objective - the objective function
% IN  : coefficients alpha and data y
% OUT : value of f(alpha,y)+valid(alpha,y)
   val = f(DzY(alpha,y))+valid(alpha,y);
end

function b = readyOrNot(values,flag)
% readyOrNot - helps to decide to break the computation
% IN  : values of objective, flag. Valid flags are: 'warn' or 'stopp'
    b = false;
    ready = all(abs(values(1:10)-values(1))<1e-6);
    if ready
        switch flag
            case 'warn'
                warning('readyOrNot:warning',...
                    ['The first ten values are nearly the same,',...
                     ' maybe we should stopp here.'])
            case 'stopp'
                warning('readyOrNot:warning',...
                    ['The first ten values are nearly the same,',...
                    ' I stopp here.'])
                b = true;
            otherwise
                MException('readyOrNot:WrongFlag',...
                    ['The flag ',char(flag),' is not valid']).throw;                              
        end
    end        
end