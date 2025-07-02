function plotMoleFractions(obj,arg,marker)
%plotMoleFractions(obj[,n[,marker]])
%   plots the mole fractions over the time or 1D space
% The optional argument can be  a cell
% array containing the species names or a vektor of indecies of the
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

% switch obj.reactor
%     case {'tibox' 'constVolumeReactor' 'commonReactor' 'constPressReactor' 'constPressNetReactor'...
%             'reactorNet' 'constVolumeNetReactor'}
%        obj.space = obj.space;  
%     case 'counterFlowFlame'
%        obj.space = obj.space; 
%     otherwise
%         warning('plotMoleFractions:unknownReactor','The reactor is unknown. Assume a 0D problem. Good luck.')
% end
        


   
if nargin==3
    % the argument can be the number of species to be plottet or a cell
    % array containing the species names or a vektor of indecies of the
    % species to be plottet...    
    indx = arg; 
    if isempty(indx)
        [m] = max(obj.moleFractions,[],2);
        [~,indx] = sort(m);
        if length(indx)>=5
            indx = indx(end-5:end);  
        elseif length(indx)>=1
            % indx = indx
        else
            % error case
            warning('plotMoleFractions:NoFractions',...
                'The moleFractions property is empty, exit without plotting.')
            return
        end
        arg = 1; % make arg numeric
    end
    if isa(arg,'numeric')
        plot(obj.space,obj.moleFractions(indx,:),marker);
        h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));        
    elseif isa(arg,'cell')&&isa(arg{1},'char')
        plot(obj.space,obj.getMoleFractionsByName(arg),marker);
        h = legend(makeNiceSpeciesNames(indx));         
    elseif isa(arg,'char')
        plot(obj.space,obj.getMoleFractionsByName(arg),marker);
        h = legend(makeNiceSpeciesNames(indx));         
    else
        warning('plotMoleFractions:WrongInput','The input must be a cell or a vektor of indices')
        %
    end
%     set(h,'Interpreter','TeX','Box', 'on','FontSize',10,'Location','Best');
elseif nargin==2 % call with species specifies
    % the argument can be the number of species to be plottet or a cell
    % array containing the species names or a vektor of indecies of the
    % species to be plottet...
    indx = arg;   
    % single moleFractions Object
    plot(obj.space,obj.moleFractions(obj.getSpeciesIndexByName(arg),:));
    h = legend(makeNiceSpeciesNames(indx)); 
    set(h,'Interpreter','TeX','Box', 'on','FontSize',20,'Location','Best');
    xlabel('t in sec','FontSize',20); 
    title('MoleFractions','FontSize',20);
    ylabel(' X','FontSize',20);
    set(gca,'FontSize',20);    
else % pure call
    [m] = max(obj.moleFractions,[],2);
    [~,indx] = sort(m);
    indx = indx(end-5:end);
    plot(obj.space,obj.moleFractions(indx,:));
    h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));
%     set(h,'Interpreter','TeX','Box', 'on','FontSize',20,'Location','Best');
end
title('MoleFractions','FontSize',20);
xlabel('T [s]','FontSize',20); 
ylabel('X','FontSize',20);
set(gca,'FontSize',20);







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

