function individualsListOut = XingR(~,individualsListIn)
% crossing the chromosomes
% works on a list of individuals

if ~isa(individualsListIn,'individual')
    error('Xing:InPutMustBeIndividual','The input must be a cell array. Each entry must be a individual object')
end

n = length(individualsListIn);
if n==1,
    error('Xing:WrongNumberElementsInList','X-ing makes only sence if there are at least TWO individuals.')
end

m = 1;  % counter for new generation
for k = 1:n-1
    for l = k+1:n
        % we call now the x-ing method from individual class
         [individualsListOut(m),individualsListOut(m+1)] =...   
             individualsListIn(k).xingR(individualsListIn(l)); 
         m = m+2;
    end
end

end

