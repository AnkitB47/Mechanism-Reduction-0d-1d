 function makeSCENGENFile(obj,chem)
 %makeSCENGENFile(obj)
 % writes the SCENGEN file for use with DRG, SEMCM,... etc.
 
 % check if the chem is correct
 if isempty(chem.simulationTime)
     error('makeSCENGENFile:chemistryObjectFailed','No simulation times defined.') 
 end
 
 % the file open stuff 
 fid = makeFileStructureForSCENGEN;
% on the screen 
%  fid = 1;
 % now we write the content
 fprintf(fid,'SIMUTIMES !START AND DURATION(SECONDS)\n');
 fprintf(fid,['  ',num2str(min(chem.simulationTime)),'  ',num2str(max(chem.simulationTime)),'\n']);
 
 fprintf(fid,'INTPARAMS ! (H0/S;RTOL;ATOL/(MOL/CM^3))\n');
 fprintf(fid,'  1E-10  1E-4  1E-12\n');
 fprintf(fid,'OUTPUTPARAMS !START(SECONDS) AND NUMBER OF OUTPUTS (+LOG.DISTR./-LIN.DIST.)\n');
 fprintf(fid,'  1E-8       -800\n');

 fprintf(fid,'PRESSURE\n');
 fprintf(fid,['  ',num2str(length(obj.pressures)),' ! NUMBER OF PRESSURES \n']);
 fprintf(fid,['  ',num2str(obj.pressures),'\n']); 
 
 fprintf(fid,'TEMPERATURE\n');
 fprintf(fid,['  ',num2str(length(obj.temperatures)),' ! NUMBER OF TEMPERATURES \n']);
 fprintf(fid,['  ',num2str(obj.temperatures),'\n']);
 

 
 
 % we check the names of the species: It should be possible that one test
 % very different mixtures instead of only variations of one set of initial
 % species
 species = getSpeciesFromInitialFractions(obj.moleFractions{1});
 moles(:,1) = getMoleFractionsFromInitialFractions(obj.moleFractions{1});
 
 % try to make a mixture matrix
 try
     for k = 2:length(obj.moleFractions)    
         moles(:,k) = getMoleFractionsFromInitialFractions(obj.moleFractions{k});
     end       
 catch ME
     fprintf(['makeSCENGENFile:DimensionError','I catch the error ',ME.identifier,'\n']);
     fprintf('Check your scenario object\n');
 end
 fprintf(fid,'FUELMIXTURE\n');
 fprintf(fid,[num2str(length(species)),'   ',num2str(length(moles(1,:))),'  !NUMBER OF INITIALIZED REACTIVE SPECIES, NUMBER OF COMPOSITIONS\n']);
 if length(species)==length(moles(:,1))
     for k = 1:length(species)
         fprintf(fid,[species{k},'  ',num2str(moles(k,:)),'\n']);
     end
 else
     error('makeSCENGENFile:DimensionError','The number of species and the number of concentrations are not the same.')
 end
 fprintf(fid,'DILUTION\n');        
 fprintf(fid,'1 ! NUMBER OF DILUTIONS\n');
 fprintf(fid,'0 !MOLE_FRACTION OF DILUENT\n');
 
 % we need ONE diluent, try on of these species...we prefer
 

 diluent = lookForDiluent(chem);
 
 fprintf(fid,'DILUENT\n');        
 fprintf(fid,'1 1 !NUMBER OF DILUENT SPECIES, NUMBER OF COMPOSITIONS\n');
 fprintf(fid,[diluent,' 0\n']);        
 
 
 fprintf(fid,'FULLMECH\n');
 fprintf(fid,'mechanism ! FULL MECHANISM IN CHEMKIN FORMAT (*.INP)\n\n');
 st = fclose(fid);
 if st~=0
     warning('makeSCENGENFile:FCloseFailed','Closing the file failed for some unknown reasons, sorry.')
 end

 
 
 
 %  ==========================local function============================
        
     function diluent = lookForDiluent(chem)
         list = getSpeciesFromMech(chem.mechanism);
         candidates = {'AR', 'NE','HE'};
         for l = 1:3,
             [b,k] = isInListI(list,candidates{l});
             if b
                 diluent = list{k};
                 break
             end
         end
     end
 
     function fid = makeFileStructureForSCENGEN
         [success,message,messageID] = mkdir('mtibox/IN/GENE'); 
         % clean the folder
         delete mtibox/IN/GENE/*;      
         if ~success
             error(messageID,message)
         end
         fid = fopen('mtibox/IN/GENE/scenario.gen','w','n','latin1');
         
         if fid>0
             fprintf('create file ~/mtibox/IN/GENE/scenario.gen.\nReady to run scengen.\n')
         else
             error('makeSCENGENFile:fopenFailed','Canot open the file scenario.gen.')
         end
     end

 end