function varargout = isInList(list,str)
% bool = isInList(list,str) returns true if the 
% string str ist found in the cell array list. Case sentitive! 
b = false;
for kk = 1:length(list)
    b = strcmp(list{kk},str);
    if b        
        break
    end
end
switch nargout
    case 1
        varargout{1} = b;
    case 2
        varargout{1} = b;
        varargout{2} = kk;
    otherwise
        varargout = [];
end
end