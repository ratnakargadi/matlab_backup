clear all;
close all;
clc;
ARGO_dir = '/projects/bobble/Data';
ARGO_file_list = [ARGO_dir filesep 'ARGOfiles_list.txt'];
file = fopen(ARGO_file_list,'r');
count = 1;
while(~feof(file))
filename = fgetl(file);
try
time = datenum(ncread(filename,'juld_location')) + datenum('1950-01-01');
catch
time = datenum(ncread(filename,'JULD_LOCATION')) + datenum('1950-01-01');
end
if(time<datenum('2016-07-05'))
p{count} = filename;
count = count + 1;
end
end
fclose(file);

OA_files_argo = '/projects/bobble/Data/OAfiles_list.txt';
file = fopen(OA_files_argo,'wt');
for i=1:length(p)
    fprintf(file,'%s',string(p{i}));
    fprintf(file,'\n');
end
fclose(file);