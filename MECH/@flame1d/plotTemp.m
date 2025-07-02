function plotTemp(obj,optional)
% plots the temperarture over time
switch nargin
    case 1
        plot(obj.space,obj.temp);
    case 2
        try
            plot(obj.space,obj.temp,optional);
        catch ME
            fprintf(['The argument of plotTemp may be not correct,',...
                ' check the result...\n The excemtions was:\n ',ME.message,'\n'])            
            plot(obj.space,obj.temp);
        end
    otherwise
end
xlabel('t [s]','FontSize',20);
ylabel('T [K]','FontSize',20); 
set(gca,'FontSize',20);
end