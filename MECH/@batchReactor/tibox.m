function obj = tibox(obj)
%obj = tibox(obj)
% Wrapper method to call tibox from matlab.
% It creates the structures needed by tibox including the scenery file,
% solves the problem by calling tibox and read the output from the *.c file
%

% tiboxpath MUST be set in the global script

tiboxpath = '~/tiboxbin/';
fprintf('\nThis is tibox called by batchReactor.tibox\n')
fprintf('Integrates the kinetic differential equations.\n')
fprintf('Authors: Tibox(R)  was written originally by T.Nagy (C) 2010, only the Matlab wrapper was added by U.Pruefert\n')
fprintf('\n')
if isempty(tiboxpath)
    error('tibox:NoToboxPath','There is no tiboxpath set.')
end
% define the reactor field
obj.reactor = 'tibox';
obj = obj.clearFractions;

% current path. For windows maschines we use this as the base folder, For
% LINUX maschines we will use the users home ~/
cpath = pwd;

% first part: preparations
% creating a LOCAL tibox folder

% 1.1 create some directories, and copy some files...
% we use xcd, which

% cd should be system independent
if isunix
%     [success,message,messageID] = xcd('~/');
    [success,message,messageID] = xcd(cpath);
    
elseif ispc
    [success,message,messageID] = xcd(cpath);
end
if ~success
    error(messageID,message)
end
% mkdir is also  system independent
[success,message,messageID] = mkdir('mtibox');
if ~success
    error(messageID,message)
end
% copy the tibox & trans programs
if isunix   
    [success,message,messageID] = copyfile([tiboxpath,'tibox'],'mtibox');
    if ~success
        error(messageID,message)
    end
    [success,message,messageID] = copyfile([tiboxpath,'trans'],'mtibox'); 
    if ~success
        error(messageID,message)
    end
elseif ispc
    [success,message,messageID] = copyfile([tiboxpath,'tibox.exe'],'mtibox');
     
    if ~success
        error(messageID,message)
    end
    [success,message,messageID] = copyfile([tiboxpath,'trans.exe'],'mtibox'); 
    if ~success
        error(messageID,message)
    end    
else
    error('tibox:UnsupportetOS','The matlab version of tibox only supports Linux/Unix and Windows systems')
end

%preparing the folders for tibox
[success,message,messageID] = xcd('mtibox');
if ~success
    error(messageID,message)
end
[success,message,messageID] = mkdir('IN');
if ~success
    error(messageID,message)
end
[success,message,messageID] = mkdir('IN/SCEN');
if ~success
    error(messageID,message)
end
[success,message,messageID] = mkdir('OUT');
if ~success
    error(messageID,message)
end
[success,message,messageID] = mkdir('OUT/CONC');
if ~success
    error(messageID,message)
end

% we need a copy of the mechanism in the tibox folder: We create as doft
% link to the mechanism file an name it mechanism.inp
% THIS MUST BE A CHEMKIN FILE!!!
% 

if isunix
    [status,result] = system(['ln -s -f ',obj.mechanism(1:end-3),'inp',' mechanism.inp']);
elseif ispc
    [status,result] = system(['copy /y ',obj.mechanism(1:end-3),'inp',' mechanism.inp']);
else
    error('tibox:UnsupportetOS','The matlab version of tibox only supports Linux/Unix and Windows systems')
end
    
if status
    error('tibox:error',result);    
else
    if isunix
        fprintf(['tibox: linking mechanism.inp ',obj.mechanism,'\n']);
    else  
        fprintf('tibox: copy mechanism.inp into the mtibox folder\n');
    end
end

% prepare the scen file for tibox
fid = fopen('IN/SCEN/scenario.scen','w','n','latin1');
 
fprintf(fid,'VERB MKDIR COPY MATRIXINP CONC KINALC JACOBI (T/F)\n');
fprintf(fid,' F T F T F F F\n');
fprintf(fid,'START(s),DURATION(s)\n');
fprintf(fid,'%19.12f %23.8f\n',obj.simulationTime);
fprintf(fid,'INTPARAMS(h0/s;rtol;atol/(mol/cm^3))\n');
fprintf(fid,'  1.000000000000000E-010  1.000000000000000E-004  1.000000000000000E-012\n');
fprintf(fid,'OUTPUT START(s), NUMBER OF OUTPUTS (+=>log,-=>lin)\n');
fprintf(fid,'  1.000000000000000E-008         -800\n');
fprintf(fid,'TEMPERATURE(K), PRESSURE(Pa)\n');
fprintf(fid,'%19.12f %23.8f\n',[obj.initialTemperature,obj.initialPressure]);
fprintf(fid,'FULL MECHANISM PREFIX(IN/INTERNAL/*. .txt)\n'); 
% mechanism w/o the file extension!
% fprintf(fid,[mechname,'\n']); % mechanism w/o the file extension
fprintf(fid,'mechanism\n');
fprintf(fid,'MECHANISM FILE(IN/REAC/*.reac or IN/MATR/*.matr)\n');
% fprintf(fid,[mechname,'\n']); 
fprintf(fid,'mechanism\n');
% this is a ~: We do not need any screen output
fprintf(fid,'SCREEN OUTPUT(list(0-10) can be in more lines, put END to the end)\n');
fprintf(fid,'END\n');
% the file output: this is the only line we can chose the output. There is
% a limitation of the number of outputs what is not in cantera...
% at least, we write all the initial species in this list
fprintf(fid,'ASCII FILE OUTPUT(list(0-15) can be in more lines, put END to the end)\n');
% fprintf(fid,selectSpeciesFromList(obj.initialMoleFractions)) 
if isempty(obj.speciesNames)
    warning('tibox:EMptySpeciesList','The List of species for output is empty, use the species from initialMoleFractions');
    listTheSpecies(fid,selectSpeciesFromList(obj.initialMoleFractions))
else
    listTheSpecies(fid,obj.speciesNames)
end
fprintf(fid,'END\n');
% the inital mole fraction, we have to normalize them,...
species = selectSpeciesFromList(obj.initialMoleFractions); % cell array with the names
moles = selectMoleFractionsFromList(obj.initialMoleFractions); % double array with the values
% etsting the data
if ~length(species)==length(moles)
    error('tibox:internalError','The number of species do not fit the number of concentrations')
end
fprintf(fid,'INITIAL MIXING RATIOS(number(1-20),species,value,metric:pct/ppm/ppb/ppt)\n');
for k = 1: length(species)
    fprintf(fid,'   %3s',species{k});    fprintf(fid,'    ');
    fprintf(fid,'   %10s',num2str(moles(k)));    fprintf(fid,'  pct            \n');
end
fclose('all');
% scenerio file is now written

% perform the the transformation of the mech file into the tibox format by
% calling trans
if isunix
    [status,result] = system('./trans mechanism');
elseif ispc
    [status,result] = system('trans mechanism');
    
else
      error('tibox:UnsupportetOS','The matlab version of tibox only supports Linux/Unix and Windows systems')
end
if status~=0
    error('tibox:error',result) 
else
    fprintf('tibox: trans runs ... ');
    fprintf(' successfully ......... 100%%\n')
    fprintf('\n')  
end


% second part: solve the problem by calling tibox. Tibox writes the result
% in a file scenario.c
if isunix
    [status,result] = system('./tibox scenario');
elseif ispc
    [status,result] = system('tibox scenario');
else
    error('tibox:UnsupportetOS','The matlab version of tibox only supports Linux/Unix and Windows systems')
end
if status~=0    
    fprintf(['tibox:system call tibox executable: error status: ',...
        num2str(status),',Message: ',result])    
else
    fprintf('tibox: tibox runs ... ');
    if strcmp(result(end-2:end),'100')
        fprintf(' successfully......... 100%%\n')
    else
        fprintf(' with problems \n')
        disp(result)
    end   
    fprintf('\n')  
end


% third part: red the *.c file and put the data in obj.moleFractions. 
 try
    data  = readTiboxC('OUT/CONC/scenario.c');
    obj.times =  data.data(1,:);
    obj.pres = data.data(2,:);
    obj.temp = data.data(3,:);
    % data  contains possibly doubled species, caused by an error in
    % tibox:
    index = checkList(obj.speciesNames,data.title);   
    obj = obj.setMoleFractions(data.data(index,:));   
    obj = obj.setMassFractions([]);
    % the ignition time is defined as argmax(dTemp/dt)
    [~,indx] = max(dot(obj.times,obj.temp)); 
    obj.ignt = obj.times(indx);
    % clear the mass fractions: tibox do not provide them    
    obj.dens = [];
catch ME
    % 
    fprintf('tibox: ')
    fprintf(ME.message);
end
% [status,result] = system('rm OUT/CONC/scenario.c');

if isunix
    [success,message,messageID] = xcd('~/');
elseif ispc
    [success,message,messageID] = xcd(cpath);
end
if ~success
    error(messageID,message)
end


%

end
% some  function, well known from e.g. spasityPattern, but with some
% adaptions
    function speciesList = selectSpeciesFromList(str) 
        % extracts species names from the (user input) inital mole fractions
        % here we use it for selecting the three body from their definiton
        % sequence what is quite similar to the cantera systax of the mole
        % fractions
        speciesList = [];
        k = 1; % zaehler fuer die position
        in = true;
        species = [];
        for i=1:length(str)        
            if strcmp(str(i),':')
                in = false; 
                speciesList{k} = species;
                species = [];
                k = k+1; % next  
            elseif isLetter(str(i))
                in = true;
            end
            if allowedChar(str(i)) && in
                species = [species,str(i)];             
            end
        end
    end

% 
    function moles = selectMoleFractionsFromList(str)
        % extracts molefractions from the (user input) inital mole
        % fractions and normalze them (if necessary)
        speciesMole =[];
        k = 1; 
        in = false;
        for i=1:length(str)   
            if strcmp(str(i),':')
                in = true;               
            end
            if isnumber(str(i)) && in
                speciesMole = [speciesMole,str(i)];             
            end
            if (strcmp(str(i),',')||i==length(str)) && in                
                moles(k) = str2double(speciesMole);
                k=k+1;
                in = false;
                speciesMole=[];
            end
        end
        moles = moles/sum(moles);
    end

    function b = isnumber(char)
        % checks if char is a number, i.e. it contains '0'...'9' and '.'
        b = ((double(char)<58 && double(char)>47))||double(char)==46;
    end

    function b = allowedChar(charSpec)
    % is it a allowed character in the definietion of a species?   
    b = (double(charSpec)<91 && double(charSpec)>64) ||...% capitals A...Z
                (double(charSpec)<123 && double(charSpec)>96) ||...lowercasses a..z
                (double(charSpec)<59 && double(charSpec)>47) ||... % numbers 0,1...9
                double(charSpec)==44||... % ,
                double(charSpec)==95||... % _
                double(charSpec)==35||... % -
                double(charSpec)==45||... % -
                double(charSpec)==43||... % +
                double(charSpec)==41||... % (
                double(charSpec)==40;     % )               
    end
    
    function b = isLetter(charSpec)
        % is it a  letter character (to define a species name, the first must
        % be a letter, 20H is NOT a well formed species name
        b = (double(charSpec)<91 && double(charSpec)>64) ||...% capitals A...Z
            (double(charSpec)<123 && double(charSpec)>96);  % lowercasses a..z
    end

    function listTheSpecies(fid,list)
    m = ceil(length(list)/8);
    r = length(list);
    n = 1;
    
    for k = 1:m,
        for l = 1:8,
            if n <= r
                fprintf(fid,[list{n},'   ']);
                n = n+1;
            end
        end
        fprintf(fid,'\n');
    end
    end
    
%     function name = fileNameFromPath(mechanismPath)
%     j = 0;
%     m = 0;
%     for k = 1:length(mechanismPath)
%         if strcmp(mechanismPath(k),'/')
%             j = j+1; 
%             slash(j)=k; 
%         elseif strcmp(mechanismPath(k),'.')
%             m=m+1;
%             dot(m)=k; 
%         end
%     end
%     name = mechanismPath(max(slash)+1:max(dot)-1);    
%     end     
    
    function list = checkList(listA,listB)
        % checks if a element of  listB is in listA and gives back the index
        m = length(listA);
        n = length(listB);
        list = 1:m;
        for k = 1:m
            for l = 1:n
                if strcmp(listA(k),listB(l))
                    list(k) = l;                   
                    break;
                end
            end
        end        
    end
    
    function dt = dot(time,value)
        delta_t = time(2:end)-time(1:end-1);
        dt = (value(:,2)-value(:,1))./delta_t(1);
        dt = [dt,(value(:,3:end)-value(:,1:end-2))*diag(1./(delta_t(1:end-1)+delta_t(2:end)))];
        dt = [dt,(value(:,end)-value(:,end-1))./delta_t(end)];
    end
    
    
   