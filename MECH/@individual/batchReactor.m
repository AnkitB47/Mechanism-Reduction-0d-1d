function batchRe = batchReactor(obj)
% batchReactor(mechanism)
if nargin == 0
    error('individual:batchReactor',...
        'Creating  batchReactor from individual Object failed.')
end

user = getenv('USER');
[success,message,~] = mkdir(['/home/',user],'tempdir') ;
if success < 1 && ~strcmp(messageID,'MATLAB:MKDIR:DirectoryExists')     
        error('ga:evaluate:mkdirFailed',['Making the temporary directory failed with ',...           
        message]);     
end

% checks, what file we have in the reference object
mechFileExtension = obj.mechanism(max(strfind(obj.mechanism,'.'))+1:end);

mechFile = ['/home/',user,'/tempdir/tempMechfile',num2str(randi(999,1,1)),'.',mechFileExtension];
            % first make the mechanism        
           
batchRe = batchReactor();
batchRe.removeSpecies(obj.mechanism,mechFile,obj.speciesNotInclude,...
                obj.names(obj.speciesNotInclude));
end   