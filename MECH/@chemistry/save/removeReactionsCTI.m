function removeReactionsCTI(mechFileIn,mechFileOut,varargin)

if strcmp(mechFileIn,mechFileOut)
    error('removeReactionsCTI:SameInputOutputFile','The input and the output files are the same, are you sure? However, it is not allowed. Bye!')
end

blocks = indexCtiFile(mechFileIn);

% Check whether any index array is empty
% If there are no reactions, three-bodies or falloffs, we set the line to 1e5.
% Then the mechanism is already a non-reaction mechanism and we will copy
% it. reactionsStartAt is now larger then the number of lines.
if isempty(blocks.reactions)
    reactions = 1e5;
else
    reactions = blocks.reactions;
end
if isempty(blocks.three_bodies)
    three_bodies = 1e5;
else
    three_bodies = blocks.three_bodies;
end
if isempty(blocks.falloffs)
    falloffs = 1e5;    
else
    falloffs = blocks.falloffs;
end

% a boolean for swithcind all <=> single reaction remove
removeSingleReactions = nargin>2;
isToRemove = [];
if removeSingleReactions
    indexSet = varargin{1};
    % create the vector of lines to be removed
    allReactions = sort([blocks.reactions blocks.three_bodies blocks.falloffs]);
    
    for k = 1:length(indexSet)
        if indexSet(k)+1>length(allReactions)
            isToRemove = [isToRemove allReactions(indexSet(k)):allReactions(indexSet(k))+6];
        else
            isToRemove = [isToRemove allReactions(indexSet(k)):allReactions(indexSet(k)+1)-1];
        end
    end
end



reactionsStartAt = min([reactions(1) ...
    three_bodies(1) falloffs(1) ]);

% tries to open the files
[fidMech,mesg] = fopen(mechFileIn,'r');
if fidMech < 0 
    fprintf(mesg)
end 
[fidNewMech,mesg] = fopen(mechFileOut,'w');
if fidNewMech < 0 
    fprintf(mesg)
end

% some comment add to the mechanism file
fprintf(fidNewMech,['# This is a reduced mechanism based on ',mechFileIn,' \n']);
fprintf(fidNewMech,'# The removing process is performed by removeReactionsCTI\n');
fprintf(fidNewMech,'# \t\t (c) U. Pruefert 2011.\n\n');

nLine = 0;

while 1 % loop over all lines
    cline = fgetl(fidMech); 
    nLine = nLine + 1;
    if ~ischar(cline)
        break
    end
    if (~ischar(cline)||nLine==reactionsStartAt-1)&&~removeSingleReactions,         
        break,
    elseif isIn(nLine,isToRemove)
        % do nothing...
    else
        fprintf(fidNewMech,[cline,'\n']);
    end
end

% we close all files
err = fclose(fidNewMech);
if err < 0
    fprintf('An error during closing the file is occured. Please check the file names etc.\n')
end
err = fclose(fidMech);
if err < 0
    fprintf('An error during closing the file is occured. Please check the file names etc.\n')
end

    function b = isIn(x,y)
        % true if x is in y
        if isempty(y)
            b = false;
            return
        end
        for kk = 1:length(y)
            b = x==y(kk);
            if b
                break
            end
        end
    end

end