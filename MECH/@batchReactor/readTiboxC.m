function   a  = readTiboxC(pathToFile)
%a  = readTiboxC(pathToFile)
% Imports a tibox output *.c-file and converts the strings to numericals.
% In the 1st. line, the header of the "concentration" table is assumed. All
% date is written in line-format, i.e. the time vector and the
% "concentrations" are the lines of the a.data field.

fid = fopen(pathToFile,'r');
titleline = fgetl(fid);
j = 1;
in = false;
stringE = [];
for k = 1:length(titleline)
    if ~strcmp(titleline(k),' ')
        stringE = [stringE,titleline(k)];
        in = true;
    elseif strcmp(titleline(k),' ')&& in
        title{j} = stringE;
        stringE = [];
        in = false;
        j = j+1;
    end
    title{j} = stringE; % the last in the list
end

a.title = title;

k = 1;
while 1
    cline = fgetl(fid);
    if ~ischar(cline)
        break
    else
        a.data(:,k) = str2num(cline); %#ok<ST2NM>
        k = k+1;
    end
end

st = fclose('all');
if st == -1
    error('readTiboxC:ErrorClosingFIle','An error occcues when closing the c-file')
end

end

