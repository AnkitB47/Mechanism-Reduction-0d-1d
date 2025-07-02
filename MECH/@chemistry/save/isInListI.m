function [b,varargout] = isInListI(list,str)
% bool = isInList(list,str) returns true if the 
% string str ist found in the cell array list. NOT case sensitive!!!
b = false;
for kk = 1:length(list)
    b = strcmpi(list{kk},str);
    if b        
        break
    end
end
switch nargout
    case {0 1}
        % b = b :-)
    case 2       
        varargout{1} = kk;
    otherwise
        error('isInListI:WrongNumberOutputs','The number of outputs can be one or two.')
end
end