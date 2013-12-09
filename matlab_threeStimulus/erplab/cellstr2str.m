% By Javier Lopez-Calderon
function st = cellstr2str(cst)

if ~iscellstr(cst)
   error('Input must be a celltring')     
end

n = length(cst);
st = '{';
for i=1:n
        st = sprintf('%s ''%s'' ', st, cst{i});
end
st = sprintf('%s}', st);