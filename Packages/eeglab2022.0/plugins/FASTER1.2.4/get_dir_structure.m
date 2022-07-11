function folders=get_dir_structure(current,base)
x=filesep;
if x=='\'
    x='\\';
end

Current=textscan(current,'%s','Delimiter',x);
Current=Current{1};
Base=textscan(base,'%s','Delimiter',x);
Base=Base{1};
[x y]=setdiff(Current,Base);
[y1]=sort(y);
folders={Current{y1}};