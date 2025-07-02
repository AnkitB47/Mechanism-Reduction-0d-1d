function plotMassFractions(obj,arg,marker)
%plotMassFractions(obj[,n[,marker]])
%   plots the mass fractions over the time or 1D space
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
%        obj.times = obj.times;  
%     case 'counterFlowFlame'
%        obj.times = obj.space; 
%     otherwise
%         warning('plotmassFractions:unknownReactor','The reactor is unknown. Assume a 0D problem. Good luck.')
% end

  


   
if nargin==3
    % the argument can be the number of species to be plottet or a cell
    % array containing the species names or a vektor of indecies of the
    % species to be plottet...    
     
    if isempty(arg)
        [m] = max(obj.massFractions,[],2);
        [~,indx] = sort(m,'descend');
        if length(indx)>=5
           indx = indx(1:6);
        elseif length(indx)>=1
            % indx = indx
        else
            % error case
            warning('plotMassFractions:NoFractions',...
                'The massFractions property is empty, exit without plotting.')
            return
        end
        arg = 1; % make arg numeric
    elseif isnumeric(arg)
         [m] = max(obj.massFractions,[],2);
         [~,indx] = sort(m,'descend');
         indx = indx(1:arg);
    end
    if isa(arg,'numeric')
        f = semilogx(obj.times,obj.massFractions(indx,:),marker);
        h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));        
    elseif isa(arg,'cell')&&isa(arg{1},'char')
        f = semilogx(obj.times,obj.getMassFractionsByName(arg),marker);
        h = legend(makeNiceSpeciesNames(arg));     
         
    elseif isa(arg,'char')
        f = semilogx(obj.times,obj.getMassFractionsByName(arg),marker);
        h = legend(makeNiceSpeciesNames(arg));         
    else
        warning('plotMassFractions:WrongInput','The input must be a cell or a vektor of indices')
        %
    end
%     set(h,'Interpreter','TeX','Box', 'on','FontSize',10,'Location','Best');
elseif nargin==2 % call with species specifies
    % the argument can be the number of species to be plottet or a cell
    % array containing the species names or a vektor of indecies of the
    % species to be plottet...
     
    if isempty(arg)
        [m] = max(obj.massFractions,[],2);
        [~,indx] = sort(m,'descend');
        if length(indx)>=5
           indx = indx(1:6);
        elseif length(indx)>=1
            % indx = indx
        else
            % error case
            warning('plotMassFractions:NoFractions',...
                'The massFractions property is empty, exit without plotting.')
            return
        end
        arg = 1; % make arg numeric
    elseif isnumeric(arg)
         [m] = max(obj.massFractions,[],2);
         [~,indx] = sort(m,'descend');
         indx = indx(1:arg);         
    end
    if isa(arg,'numeric')
        f = semilogx(obj.times,obj.massFractions(indx,:));
        h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));        
    elseif isa(arg,'cell')&&isa(arg{1},'char')
        f = semilogx(obj.times,obj.getMassFractionsByName(arg));
        h = legend(makeNiceSpeciesNames(arg));         
    elseif isa(arg,'char')
        f = semilogx(obj.times,obj.getMassFractionsByName(arg));
        h = legend(makeNiceSpeciesNames(arg));         
    else
        warning('plotMassFractions:WrongInput','The input must be a cell or a vektor of indices')
        %
    end
%     set(h,'Interpreter','TeX','Box', 'off','FontSize',18,'Location','Best');
else % pure call
    [m] = max(obj.massFractions,[],2);
    [~,indx] = sort(m,'descend');
    indx = indx(1:min(6,size(obj.massFractions,1))); 
    f = semilogx(obj.times,obj.massFractions(indx,:));
    h = legend(makeNiceSpeciesNames(obj.speciesNames(indx)));
    
end
title('MassFractions','FontSize',18);
xlabel('t [s]','FontSize',18); 
ylabel('Y','FontSize',18);
set(gca,'FontSize',18);
set(f,'LineWidth',2);
set(h,'Interpreter','TeX','Box', 'on','FontSize',18,'Location','Best','EdgeColor','w');





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

