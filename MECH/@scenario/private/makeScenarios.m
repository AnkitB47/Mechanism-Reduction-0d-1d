function makeScenarios(obj,chem) 
%makeScenarios(obj,chem) 
% Writes a set of tibox scenario files
% from the data includet in the scenario object
% It is of course a privat method and will be called exclusevely by
% makeScenarioFiles.

% create the folder mtibox/IN/SCEN

[success,message,messageID] = mkdir('mtibox/IN/SCEN/'); 
% clean the folder

if ~success
    error(messageID,message)
end
delete mtibox/IN/SCEN/*;
% link the base mechanism
% pathToMech = chem.mechanism;  
% [status,result] = system(['ln -s ',chem.mechanism,' ',...
%     'mtibox/',pathToMech(max(strfind(pathToMech,'/'))+1:end)]);
% if status~=0 
%     warning('makeScenarios:systemCall',result);
% end

nOfScenarios = length(obj.temperatures)*...
    length(obj.pressures)*length(obj.moleFractions);

% three loops over the scenarios
msg(['Try to write ',num2str(nOfScenarios),' scenarios']);
for k = 1:length(obj.temperatures)
    for kk = 1:length(obj.pressures)
        for kkk = 1:length(obj.moleFractions);
            filename = ['mtibox/IN/SCEN/scenario_temp_',...
                num2str(obj.temperatures(k)),'_press_',...
                num2str(obj.pressures(kk)),'_conc_',num2str(kkk),'.scen'];
            fid = fopen(filename,'w','n','latin1');
            if fid>0
                msg(['create file ', filename,'\n']);
            else
                error('makeSCENGENFile:fopenFailed','Canot open the file scenario.gen.')
            end
            writeScenarioFiles(fid,chem.simulationTime,obj.temperatures(k),...
                obj.pressures(kk),obj.moleFractions{kkk});
            status = fclose(fid);
            if status~=0
                msg('fclose failed');
            end
        end
    end
end
            





    function writeScenarioFiles(fid,times,temp,press,moleFracs)
        fprintf(fid,'VERB MKDIR COPY MATRIXINP CONC KINALC JACOBI (T/F)\n');
        fprintf(fid,' F T F T F F F\n');
        fprintf(fid,'START(s),DURATION(s)\n');
        fprintf(fid,'%19.12f %23.8f\n',times);
        fprintf(fid,'INTPARAMS(h0/s;rtol;atol/(mol/cm^3))\n');
        fprintf(fid,'  1.000000000000000E-010  1.000000000000000E-004  1.000000000000000E-012\n');
        fprintf(fid,'OUTPUT START(s), NUMBER OF OUTPUTS (+=>log,-=>lin)\n');
        fprintf(fid,'  1.000000000000000E-008         -800\n');
        fprintf(fid,'TEMPERATURE(K), PRESSURE(Pa)\n');
        fprintf(fid,'%19.12f %23.8f\n',[temp,press]);
        fprintf(fid,'FULL MECHANISM PREFIX(IN/INTERNAL/*. .txt)\n'); 
        % mechanism w/o the file extension!
        % fprintf(fid,[mechname,'\n']); % mechanism w/o the file extension
        fprintf(fid,'mechanism\n');
        fprintf(fid,'MECHANISM FILE(IN/REAC/*.reac or IN/MATR/*.matr)\n');
        fprintf(fid,'mechanism\n');
        % this is a ~: We do not need any screen output
        fprintf(fid,'SCREEN OUTPUT(list(0-10) can be in more lines, put END to the end)\n');
        fprintf(fid,'END\n');
        
        % The file output: this is used for computing the error.
        % We write all the initial species in this list. 
        fprintf(fid,'ASCII FILE OUTPUT(list(0-15) can be in more lines, put END to the end)\n');
        species = getSpeciesFromInitialFractions(moleFracs);
        for l = 1: length(species)
            fprintf(fid,'  %4s ',species{l});fprintf(fid,'\n');
        end
        fprintf(fid,'END\n');
        
        % Nowe we write the initial mole fractions
        moleFractions = getMoleFractionsFromInitialFractions(moleFracs);
       
        % the inital mole fraction, we have to normalize them,...
        % teststing the data
        if ~length(species)==length(moleFractions)
            error('tibox:internalError','The number of species do not fit the number of concentrations')
        end
        fprintf(fid,'INITIAL MIXING RATIOS(number(1-20),species,value,metric:pct/ppm/ppb/ppt)\n');
        for l = 1: length(species)
            fprintf(fid,'   %4s',species{l});    fprintf(fid,'    ');
            fprintf(fid,'   %10s',num2str(moleFractions(l)));    fprintf(fid,'  pct            \n');
        end
        fprintf(fid,'\n');
    end




    function msg(text)
        fprintf('makeScenarios:\n\t');
        fprintf(text);fprintf('\n');
    end



end

