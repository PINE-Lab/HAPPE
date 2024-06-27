fid = fopen('test_str2array.csv','r');	%% corrected directory 
if fid<0, return; end; 
s   = fread(fid,[1,inf],'uint8=>char');
fclose(fid); 
s(s==10)=[];
[n,v,c]=str2array(s,[';',char(9)],char([10,13]))


