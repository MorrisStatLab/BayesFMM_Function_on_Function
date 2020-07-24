function Zcat=cat_cell(Z)
a=length(Z);
Zcat=[];
for i=1:a
    Zcat=[Zcat,Z{i}]; %#ok<AGROW>
end