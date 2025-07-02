function plotMoleFractions(obj,arg,marker)
%plotMoleFractions(obj[,n[,marker]])
%   plots the mole fractions over the time or 1D space
% The optional argument the species names or a vektor of indecies of the
% species to be plottet.. If no argument is given, the six "most active" species
% are selected for plotting.
% The optional argument marker can be every charakter from the following
%  b     blue          .     point              -     solid
%  g     green         o     circle             :     dotted
%  r     red           x     x-mark             -.    dashdot 
%  c     cyan          +     plus               --    dashed   
%  m     magenta       *     star             
%  y     yellow        s     square
%  k     black         d     diamond
%  w     white         v     triangle (down)
%                      ^     triangle (up)
%                      <     triangle (left)
%                      >     triangle (right)
%                      p     pentagram
%                      h     hexagramsuitable 

  


   
if nargin==3
    % the argument can be the number of species to be plottet or a cell
    % array containing the species names or a vektor of indecies of the
    % species to be plottet...    
     
    if isempty(arg)
        [m] = max(obj.moleFractions,[],2);
        [~,indx] = sort(m,'descend');
        if length(indx)>=5
           indx = indx(1:6);
        elseif length(indx)>=1
            % indx = indx
        else
            % error case
            warning('plotMoleFractions:NoFractions',...
                'The moleFractions property is empty, exit without plotting.')
            return
        end
        arg = 1; % make arg numeric
    elseif isnumeric(arg)
         [m] = max(obj.moleFractions,[],2);
         [~,indx] = sort(m,'descend');
         indx = indx(1:arg);
    end
    if isa(arg,'numeric')
        plot(obj.times,obj.moleFractions(indx,:),marker);
        h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));        
    elseif isa(arg,'cell')&&isa(arg{1},'char')
        plot(obj.times,obj.getMoleFractionsByName(arg),marker);
        h = legend(makeNiceSpeciesNames(arg));     
         
    elseif isa(arg,'char')
        plot(obj.times,obj.getMoleFractionsByName(arg),marker);
        h = legend(makeNiceSpeciesNames(arg));         
    else
        warning('plotMoleFractions:WrongInput','The input must be a cell or a vektor of indices')
        %
    end
    set(h,'Interpreter','TeX','Box', 'on','FontSize',10,'Location','Best');
elseif nargin==2 % call with species specifies
    % the argument can be the number of species to be plottet or a cell
    % array containing the species names or a vektor of indecies of the
    % species to be plottet...
     
    if isempty(arg)
        [m] = max(obj.moleFractions,[],2);
        [~,indx] = sort(m,'descend');
        if length(indx)>=5
           indx = indx(1:6);
        elseif length(indx)>=1
            % indx = indx
        else
            % error case
            warning('plotMoleFractions:NoFractions',...
                'The moleFractions property is empty, exit without plotting.')
            return
        end
        arg = 1; % make arg numeric
    elseif isnumeric(arg)
         [m] = max(obj.moleFractions,[],2);
         [~,indx] = sort(m,'descend');
         indx = indx(1:arg);         
    end
    if isa(arg,'numeric')
        plot(obj.times,obj.moleFractions(indx,:));
        h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));        
    elseif isa(arg,'cell')&&isa(arg{1},'char')
        plot(obj.times,obj.getMoleFractionsByName(arg));
        h = legend(makeNiceSpeciesNames(arg));         
    elseif isa(arg,'char')
        plot(obj.times,obj.getMoleFractionsByName(arg));
        h = legend(makeNiceSpeciesNames(arg));         
    else
        warning('plotMoleFractions:WrongInput','The input must be a cell or a vektor of indices')
        %
    end
    set(h,'Interpreter','TeX','Box', 'on','FontSize',10,'Location','Best');
else % pure call
    [m] = max(obj.moleFractions,[],2);
    [~,indx] = sort(m,'descend');
    indx = indx(1:6); 
    plot(obj.times,obj.moleFractions(indx,:));
    h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));
    set(h,'Interpreter','TeX','Box', 'on','FontSize',10,'Location','NorthEastOutside');
end
title('MoleFractions','FontSize',14);
xlabel('T [s]','FontSize',14); 
ylabel('X','FontSize',14);
set(gca,'FontSize',14);







    function str = makeNiceSpeciesNames(str)
        if isa(str,'cell')
            for k = 1:length(str)
                            str{k} = makeSubscriptSpecies(str{k});
            end
        elseif isa(str,'char')
            str = makeSubscriptSpecies(str);
        else
                % hmmmm...?
        end
        function str = makeSubscriptSpecies(str)
            strd = [];
            l = 1;
            while l <= length(str)                
                if double(str(l))>=48 && double(str(l))<=57 
                    strd = [strd,'_{'];
                    while l <= length(str)                       
                        if double(str(l))>=48 && double(str(l))<=57 && l< length(str)
                            strd = [strd, str(l)];
                            l = l+1;  
                        elseif double(str(l))>=48 && double(str(l))<=57 && l == length(str)
                            strd = [strd,str(l),'}'];
                            l = l+1;                            
                            break
                        else
                            strd = [strd,'}'];                             
                            break
                        end
                    end
                else
                   strd = [strd,str(l)];
                   l = l+1;
                end
                
            end
            str = upper(strd);
        end
    end

end

