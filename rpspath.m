%%
% put all folders into the default path
subfolders={'assem', 'basis','ex_a','geops','gmesh','solve','utility'};
if ispc
    sc='\';
else
    sc='/';
end
 for k=1:length(subfolders)
        mypaths=strcat(pwd,sc,subfolders{k});
        path(path,mypaths)
 end 
 clear all
   
