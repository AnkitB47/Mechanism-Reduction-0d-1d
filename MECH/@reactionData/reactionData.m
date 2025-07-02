classdef reactionData
    % class to store the reaction related data   
    % reactionData(XML-FilE)
    % reactionData(XMLTREE) 
    % (c) 2011 by Uwe Pruefert vor VIRTUHCON
    properties(SetAccess=private)   
        A;
        b;
        E;
        rtype;
        efficiencies;
        units;
        
    end
    
    methods
        function obj = reactionData(varargin)
            % class to store the reaction related data   
            % reactionData(XML-FilE)
            % reactionData(XMLTREE) 
            % (c) 2011 by Uwe Pruefert vor VIRTUHCON
            % Fields are 
            % rd.A              frequence factor
            % rd.b              pre-exponetial temperature exponent
            % rd.E              activation energy
            % rd.efficiencies   stores the effciencies for tree body and
            %                   fall-off reactions
            % rd.units          gives the unit of E, e.g. cal/mol or kcal/mol
            % (depends on the mechanism)
            
    
            switch nargin
                case 0
                case 1
                    if  isstruct(varargin{1})
                        xmlTree = varargin{1};
                        if isempty(xmlTree(1).speciesData)
                            error('reactionData:NoXMLTree','The xml tree input is incorrect.')
                        end
                    elseif ischar(varargin{1})
                        try 
                            load(varargin{1});
                        catch ME
                            fprintf([ME.message,'\n']);
                            try
                                if ~exist( 'xmlTree','var')
                                    xmlTree = xml2struct(varargin{1});
                                end
                            catch ME
                                error(ME.identifier,ME.message)
                            end
                        end
                    else
                        error('reactionData:WronInputFormat','The input must be a na,e of an xml file or an xml tree');
                    end
                    
                otherwise
                    error('reactionData:WrongNumberArguments','The number of input arguments is wrong, try help reactionData');                    
            end
            
            % initialize the efficiencies matrix...
            
            obj = arrhenius(obj,xmlTree);
%             obj.A = koeffField(1,:);
%             obj.b = koeffField(2,:);
%             obj.E = koeffField(3,:);
            
            
            
            
            
            obj.units{1} = xmlTree(1).reactionData(1).reaction.rateCoeff(1).Arrhenius.E.units;
            for k2 = 2:length(xmlTree(1).reactionData)
                obj.units{2} = xmlTree(1).reactionData(k2).reaction.rateCoeff(1).Arrhenius.E.units;
                if ~strcmp(obj.units{1},obj.units{2})                    
                    error('reactionData:nonUniqueUnits',...
                        'This mechanism contains different units for its coeffiicents. This feature is not jet implemented. Sorry!')
                end
                obj.units{1} = obj.units{2};
            end
            unit = obj.units{1};
            obj.units = [];
            obj.units = unit;
                
            function obj = arrhenius(obj,xmlTree)
                % helper function to separate the coefficients to a matrix and
                % a cell array. For use with arrheniusCoeff and
                % makeArrheniusFunction
                
                % initialize the efficiencie Matrix, Only three body and
                % fall-off reaction will put contributuions to this matrix
                 
               effindex = 0;
               falloffindex = 0;     
                
                
                
                
                for k = 1:length(xmlTree(1).reactionData)
                    fd = fields(xmlTree(1).reactionData(k).reaction);
                    for kk = 1:length(fd)
                        hasReactionType = false;
                        if strcmp(fd(kk),'type')                         
                            hasReactionType = true;
                            break
                        end                        
                    end
                    if hasReactionType
                        reactionType = xmlTree(1).reactionData(k).reaction.type;                
                        % it could be a fall off or a three body                    
                        switch reactionType
                            case 'falloff'
                                % hmmm...
                                %fprintf('falloff reaction.\n')
                                effindex = effindex + 1;
                                falloffindex = falloffindex + 1;
                                % fall off reaction has two sets of
                                % arrhenius coefficients
                                obj.A(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff(1).Arrhenius.A);
                                obj.b(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff(1).Arrhenius.b);...
                                obj.E(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff(1).Arrhenius.E.E); 
                                obj.A(2,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff(2).Arrhenius.A);
                                obj.b(2,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff(2).Arrhenius.b);...
                                obj.E(2,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff(2).Arrhenius.E.E); 
                                efficienciesLine = ...
                                    xmlTree(1).reactionData(k).reaction.rateCoeff(1).efficiencies.efficiencies;
                              
                                % read the species names
                                efficienciesSpecies = getSpeciesFromInitialFractions(efficienciesLine);
                                efficienciesPerSpecies = ...
                                    getMoleFractionsFromInitialFractions(efficienciesLine,'normalize','off');
                                
                                indexEfficencieSpecies = [];
                                for k3 = 1:length(xmlTree(1).speciesData)
                                    for l= 1:length(efficienciesSpecies)
                                        if strcmp(xmlTree(1).speciesData(k3).species.name,efficienciesSpecies{l});
                                            indexEfficencieSpecies = [indexEfficencieSpecies,k3];
                                        end
                                    end
                                end
                                
                                obj.efficiencies.matrix(effindex,:) = ones(1,length(xmlTree(1).speciesData));
                                obj.efficiencies.matrix(effindex,indexEfficencieSpecies) = efficienciesPerSpecies;
                                obj.efficiencies.reactionNumber(effindex) = k;
                                obj.rtype{k}  = 'fall-off'; 
                                
                                obj.efficiencies.fallofftype{falloffindex} = xmlTree(1).reactionData(k).reaction.rateCoeff(1).falloff.type;
                                % her we have to check if the field is
                                % existent!!!
                                switch xmlTree(1).reactionData(k).reaction.rateCoeff(1).falloff.type
                                    case 'Troe'
                                        obj.efficiencies.falloff{falloffindex} = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff(1).falloff.falloff);
                                    case 'Lindemann'
                                        % nothing to do?
                                    otherwise
                                end
                                
                            case 'threeBody'
                                % three body reaction blocks provide
                                % rate_coefficients and efficiencies. The
                                % EFFICIENCIES entry is used much later to
                                % prepare the reaction line when making the
                                % w-dot function...
                               % fprintf('three body reaction.\n')
                                effindex = effindex +1;
                                obj.A(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff.Arrhenius.A);
                                obj.b(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff.Arrhenius.b);
                                obj.E(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff.Arrhenius.E.E); 
                                efficienciesLine = ...
                                    xmlTree(1).reactionData(k).reaction.rateCoeff(1).efficiencies.efficiencies;
                                 
                                % read the species names
                                efficienciesSpecies = getSpeciesFromInitialFractions(efficienciesLine);
                                % read the numbers, for efficiencies, use
                                % the option 'normalize' = 'off'
                              
                                efficienciesPerSpecies = ...
                                    getMoleFractionsFromInitialFractions(efficienciesLine,'normalize','off');
                                
                                indexEfficencieSpecies = [];

                               
                                for k3 = 1:length(xmlTree(1).speciesData)
                                    for l = 1:length(efficienciesSpecies)                                        
                                        if strcmp(xmlTree(1).speciesData(k3).species.name,efficienciesSpecies{l});  
%                                             xmlTree.speciesData(k3).species.name
%                                             efficienciesSpecies{l}
                                            indexEfficencieSpecies = [indexEfficencieSpecies,k3];
                                        end
                                    end
                                end
                                obj.efficiencies.matrix(effindex,:) = ones(1,length(xmlTree(1).speciesData));
 
                                obj.efficiencies.matrix(effindex,indexEfficencieSpecies) = efficienciesPerSpecies;
                                obj.efficiencies.reactionNumber(effindex) = k;
                                
%                                 obj.efficiencies.fallofftype{effindex}  = '-'; 
                                % tbr = Three Body Reaction
                                obj.rtype{k} = 'treebody' ;
                                
                                % ohhhh...                              
                            otherwise
                               % So what? WHAT FOR A TYPE YOU ARE?...
                               % FUCKHEAD.
                               error('reactionData:unknownType','Unkown reaction type.')
                        end
                    else
                       % no type means "common reaction, here arrhenius type reaction"
                       obj.A(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff.Arrhenius.A); %#ok<*ST2NM>
                       obj.b(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff.Arrhenius.b);
                       obj.E(1,k) = str2num(xmlTree(1).reactionData(k).reaction.rateCoeff.Arrhenius.E.E); 
                       obj.rtype{k} = 'arrhenius';
                    end                  
                end            
            end
        end
        
        
        function coeff = arrheniusCoeff(obj,temp,x,nrect)
            % k = arrheniusCoeff(T,X)
            % k = arrheniusCoeff(T,X,nreaction)
            % Computes the Arrhenius coefficient for all or for selected
            % reactions. For fall-off reactions, the concentrations are
            % neccessary.
            
            % **** Here we have to check the unit of E!!!           **** 
            % **** By computing the gas konstant in cal/mol 1/K we  ****
            % **** can now use cal/mol                              ****
            % the obj.units gives us a hint to use which unit for the gas
            % constant.
            switch obj.units
                case 'cal/mol'
                    R = 1.985877519824210;
                case 'kcal/mol'
                    R =1.985877519824210*1.0e-3;
                case 'J/mol'
                    R = 8.314472;
                case 'kJ/mol'
                    R = 8.314472*1.0e-3;
                otherwise
                    error('arrheniusCoeff:WrongUnit',['The unit ',obj.units,'is unkown. Possible  units are cal/mol, J/mol kJ/mol'])
            end
            
            switch nargin
                case 3
                    % all reactions
                    coeff = obj.A.*temp.^obj.b.*exp(-obj.E./(R*temp));
                case 4
                    coeff = obj.A(nrect)*temp.^obj.b(nrect)*exp(-obj.E(nrect)/(R*temp));
                otherwise
                    error('arrheniusCoeff:NotEnoughInputArguments','Not enough input arguments, at least call arrheniusCoeff(T,X)')
            end
            % identify fall-off reactions: There are a second non-trivial
            % row in the matrix.
            indx = find(obj.A(2,:)~=0);
            for k = 1:length(indx)
                % in this loop we correct the computation of the arrhenius
                % coeffcient for the case of fall-off reactions
                % what is x? Concentrations, Mole Fractions, ...?          
                
                M = obj.efficiencies.matrix(k,:)*x; 
                
                k_inf = coeff(1,indx(k));
                k_0 = coeff(2,indx(k));
                pr = k_0/k_inf*M;
                switch  obj.efficiencies.fallofftype{k}
                    case 'Lindemann'
                        % F = 1
                        Kf = k_inf*pr/(1+pr);
                        coeff(1,indx(k)) = Kf;
                        coeff(2,indx(k)) = 0;
                    case 'Troe'                         
                        % experimentell
                        % obj.efficiencies.falloff(1) = A
                        % obj.efficiencies.falloff(2) = T3
                        % obj.efficiencies.falloff(3) = T1
                        % obj.efficiencies.falloff(4) = T2
                        Fcent = (1-obj.efficiencies.falloff{k}(1))*exp(-temp/obj.efficiencies.falloff{k}(2))...
                            +obj.efficiencies.falloff{k}(1)*exp(-temp/obj.efficiencies.falloff{k}(3))...
                            +exp(-obj.efficiencies.falloff{k}(4)/temp);
                        C = -0.4-0.67*log10(Fcent); 
                        N = 0.75-1.27*log10(Fcent);
                        f1 = (log10(pr)+C)/(N-0.14*(log10(pr)+C));
                        log10F = log10(Fcent)/(1/(1+f1*f1));
                        F = 10^log10F;
                        Kf = k_inf*pr/(1+pr)*F;
                        coeff(1,indx(k)) = Kf;
                        coeff(2,indx(k)) = 0;
                    otherwise
                        % not  implemented types...
                        % do nothing 
%                        
                end 
                 
            end
            if ~isempty(find(coeff(2,:)~=0,1))
                error('arrheniusCoeff:ComputationFailed','One of the fall-off coefficients are not computed properly.')
            else
                coeff(2,:) = [];
            end
        end
        
        function makeArrheniusFunction(obj)
            % writes the function arrhenius(t) into a Matlab m-file
            
            fid = fopen('arrhenius.m','w');
            fprintf(fid,'function k = arrhenius(t)\n');
            fprintf(fid,'%%This is an auto generated function for the computation \n');
            fprintf(fid,'%% of the Arrhenius coefficients.\n');
            fprintf(fid,'%%\n');
            fprintf(fid,'%%Creator: reactionData.makeArrheniusFunction.\n');
            fprintf(fid,'%%(c) 2011 Uwe Pruefert for VIRTUHCON\n');
            fprintf(fid,'%%\n');
            fprintf(fid,'%%Usuage:\n%%\tk = arrhenius(t)\tt can be a double or a vector of double\n');
            fprintf(fid,'%%\n');
            fprintf(fid,'%%The result is\n%%\t (i) a row vector of dimension #reactions\n');
            fprintf(fid,'%% or\n%%\t(ii) an array of dimension length(t) x #reactions\n');
            % check the unit of E and fit the unit of the Gaskostant
            switch obj.units
                case 'cal/mol'
                    fprintf(fid,'R = 1.985877519824210; %% cal/mol * 1/K\n\n');
                case 'kcal/mol'
                    fprintf(fid,'R = 1.985877519824210e-3;%% kcal/mol * 1/K\n\n');
                case 'J/mol'
                    fprintf(fid,'R = 8.314472; %% J/mol * 1/K\n\n');
                case 'kJ/mol'
                    fprintf(fid,'R = 8.314472*1.0e-3; %% kJ/mol * 1/K\n\n');
                otherwise
                    error('arrheniusCoeff:WrongUnit',['The unit ',obj.units,'is unkown. Possible  units are cal/mol, J/mol kJ/mol'])
            end
            
            % define the coefficents. We ignore the second row in the
            % Matrizes A,b, and E, becaus they are only relavant in the
            % (separate handled) case of fall-off reactions.
            fprintf(fid,'A = [');
            for k = 1:length(obj.A)
                if obj.A(2,k)==0
                    fprintf(fid,[num2str(obj.A(1,k)),' ']);
                else
                    fprintf(fid,num2str('0 ')); 
                end
                if mod(k,8)==0
                    fprintf(fid,'...\n\t');                    
                end
            end
            fprintf(fid,'];\n');
            fprintf(fid,'\n');
            
            fprintf(fid,'b = [');
            for k = 1:length(obj.b) 
                if  obj.A(2,k)==0
                    fprintf(fid,[num2str(obj.b(1,k)),' ']);
                else
                    fprintf(fid,num2str('0 '));
                end
                if mod(k,16)==0
                    fprintf(fid,'...\n\t');                    
                end
            end
            fprintf(fid,'];\n');
            fprintf(fid,'\n');
             
            fprintf(fid,'E = [');
            for k = 1:length(obj.E) 
                if obj.A(2,k)==0 
                    fprintf(fid,[num2str(obj.E(1,k)),' ']);
                else
                    fprintf(fid,num2str('0 '));
                end
                if mod(k,8)==0
                    fprintf(fid,'...\n\t');                    
                end
            end
            
            fprintf(fid,'];\n');
            fprintf(fid,'\n');
            
            
            % print the *code* 
            % we loop over the temperature, relevant for the case that t is 
            % not a single double but a vector of temperatures
            fprintf(fid,'for l = 1:length(t)\n');
            fprintf(fid,'\tk(l,:) =  A.*t(l).^b.*exp(-E./(R*t(l)));\n'); 
            fprintf(fid,'end\n');
            
            
            % korrekt the arrhenius koefficient funtion for the fall-off
            % reactions: We loop over all fall off reations and write a
            % command like k(indx) = 
            for k = 1:lenght(obj.A(1,:))
                if strcmp(obj.rtype,'fall-off')
                    % write down the insructions to correct the k,
                    % we try to pre compute all factor as far as possible 
                else
                    % no action needed
                end
                
                
                
            end
        
            
            
            
            
            fprintf(fid,'\n');            
            fprintf(fid,'end\n');
            fclose(fid);
        end
    end
    
end

