function obj = timeScale(obj,type)
% obj = obj.timeScale(arg)
% Computes the time scale of type arg 
% arg can be pvts, gpvts, spts ijts and evts.
switch type
    case 'spts'
        obj = obj.spts();
    case 'ijts'
        obj = obj.ijts();
    case 'evts'
        obj = obj.evts();
    case 'pvts'
        obj = obj.pvts();
    case 'gpvts'
        obj = obj.gpvts();
    otherwise
        ME = MException('batchReactor:timeScale:unknownArgument',...
            ['The argument ',char(type),' is not valid.']);
        throw(ME);        
end
end

