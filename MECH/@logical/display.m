function display(obj)
fprintf([inputname(1),' =\n'])
    for k = 1:size(obj,1)
        fprintf('\t');
        for l = 1:size(obj,2)
            if obj(k,l)
                fprintf(' true ')
            else
                fprintf( 'false ')
            end
        end
        fprintf('\n');
    end
end