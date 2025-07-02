function obj = evaluate(obj,varargin)
% method for class ga to evaluate the functional obj.objective

% take the  genome an create a mechanism
% if obj.aimsToReduce 'species', we write a new mechfile
% if obj.aimsToReduce 'reactions' we manipulate the gas object by setting
% multiplier to obj.currentGeneration.reactions...

 
if ~isa(obj.currentGeneration,'individual')
    ME = MException('evaluate:wrongClasss',...
        'obj.currentGeneration must be of class ''individual''');
    ME.throw
end

if ~strcmp(obj.reference.mechanism,obj.genotype.mechanism)
    ME = MException('ga:evalutate:MechNotEqual',...
        ['The mechanisms in obj.genotype.baseMechanism and',...
        ' obj.reference.mechanism must be the same.']);
    ME.throw
end

user = getenv('USER');
[success,message,~] = mkdir(['/home/',user],'tempdir') ;
if success < 1 && ~strcmp(messageID,'MATLAB:MKDIR:DirectoryExists')     
        ME = MException('ga:evaluate:mkdirFailed',...
            ['Making the temporary directory failed with ',...           
        message]);  
    ME.throw
end

% checks, what file we have in the reference object
mechFileExtension = obj.genotype.mechanism(max(strfind(obj.genotype.mechanism,'.'))+1:end);

[fid,message] = fopen(['/home/',user,'/tempdir/results.txt'],'a');
if fid < 1
    ME = MException('ga:evaluate:fopenFailed',...
        ['Something wrong with fopen. The message was: ',message]);
    ME.throw
end
fprintf(fid,'mech                  |     f      |     nSpecies\n');
st = fclose(fid);
if st<0
    warning('ga:evaluate:FcloseFailed','Closing the data file failed');
end
tempObj = obj.reference.copy; % copy instead assigment
switch obj.aimsToReduce
    case 'species'
        for k = 1:length(obj.currentGeneration)
            mechFile = ['/home/',user,'/tempdir/mechfile',num2str(randi(1e8)),'.',mechFileExtension];
            % first make the mechanism        
            obj.reference.removeSpecies(mechFile,obj.currentGeneration(k).speciesNotInclude);
            % we check of what class the reference object is,
            % because we have to set the transport model for flames
            
            switch class(tempObj )
                case 'batchReactor'
                     tempObj.setMechanism(mechFile);
                case 'flame1d'
                     tempObj.setMechanism(mechFile,tempObj.transport);
                case 'recirculationFlame'
                     tempObj.setMechanism(mechFile,tempObj.flame.transport); 
                otherwise
                    ME = MException('ga:evaluate:unknownObject',...
                        ['The class ',class(tempObj ),...
                        ' is not known by ga.evaluate.']);
                    ME.throw
            end
            try
                % now solve and eval in one block
                tempObj.solve;    % new handle syntax                 
                obj.currentGeneration(k).value = feval(obj.objective,obj.reference,tempObj);            
            catch ME
                warning(['An error occures:',...
                    ME.message]);
                % if we catch an error, value = inf...
                obj.currentGeneration(k).value = Inf;
            end
%             obj.currentGeneration(k).value = feval(obj.objective,obj.reference,tempObj);
            
            [fid,message] = fopen(['/home/',user,'/tempdir/results.txt'],'a');
            if fid < 1
                ME = MException('ga:evaluate:fopenFailed',...
                    ['Something wrong with fopen. The message was: ',...
                    message]);
                ME.throw
            end
            fprintf(fid,[mechFile,'  |  ',num2str(obj.currentGeneration(k).value)]);
            fprintf(fid,['  |  ',num2str(tempObj.nSpecies)]);
            fprintf(fid,'\n');
            st = fclose(fid);
            if st<0
                warning('ga:evaluate:fcloseFailed','Closing the data file failed');
            end
        end
    case 'reactions'
        % not jet implemented, but comming soon 
        % The evaluation in "reactions" case works completely different.
        % We use the multiplier property of the kinetics object to switch
        % reactions on/off. The mechanism as file will not changed. 
        
        % first we try to solve the reference object
        try
             tempObj.solve; % handle syntax
        catch ME
            throwAsCaller(ME);
        end
        % now we evaluate the copy of the rerecence object
        
        for k = 1:length(obj.currentGeneration)
            tempObj.setReactions(...
                obj.currentGeneration(k).reactions);
            tempObj.solve;
            obj.currentGeneration(k).value = ...
                feval(obj.objective,obj.reference,tempObj);       
            [fid,message] = fopen(['/home/',user,'/tempdir/results.txt'],'a');
            if fid < 1
                ME = MException('ga:evaluate:fopenFailed',...
                    ['Something wrong with fopen. The message was: ',message]);
                ME.throw
            end
            fprintf(fid,num2str(obj.currentGeneration(k).value));
            fprintf(fid,['  |  ',...
                num2str(sum(obj.currentGeneration(k).reactions))]);
            fprintf(fid,'\n');
            st = fclose(fid);
            if st<0
                warning('ga:evaluate:fcloseFailed',...
                    'Closing the data file failed');
            end
        end
    otherwise        
        ME = MException('ga:evaluate:WrongTarget',...
            ['The target ',char(obj.aimsToReduce), ' is not valid']);
        ME.throw
end
end

