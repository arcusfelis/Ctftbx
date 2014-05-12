function tftbinst()

% INSTALLATION program for the ANSI C TIME-FREQUENCY TOOLBOX
% Make sure you have run mex -setup, then
% Run 'install' and follow the instructions



close all
answer=questdlg('C language Time-Frequency toolbox - Installation program. Have you read the instructions before proceeding to compilation?','Continue?');
syst_file=fullfile('src','system.h');
if strcmp(answer,'Yes'),
  choice_sys=menu('Choose your system','Windows','Unix/Linux','MacOs');
  switch choice_sys,
   case 1;3,
    system = 'windows';
    % edit the conditional compilation file
    delete(syst_file);
    fid=fopen(syst_file,'w');
    fprintf(fid,'#define SYSTEME  WINDOWS\n');
    fclose(fid);
    
   case 2,
    system = 'unix';
    % edit the conditional compilation file
    delete(syst_file);
    fid=fopen(syst_file,'w');
    fprintf(fid,'#define SYSTEME  UNIX\m\n');
    fclose(fid);
  end
  path_dest=fullfile(cd,'Ctftb');
  mkdir('Ctftb')
  
  msg=['The programs will now be compiled and the binary files' ...
       ' will be moved to '];
  msg(end+1:end+length(path_dest))=path_dest;
  answer=questdlg(msg,'Continue ?');
  if strcmp(answer,'Yes'),
    % compiles the C programs and moves them to
    % the selected directory
    tftbcomp(path_dest);
    addpath(path_dest);
    msgbox('Compilation completed !');
    warndlg('You now need to copy the directory Ctftb at its final location, and add it to your path',...
	    'One more thing to do !');
  end
else
   msgbox('Please read the installation instructions');
end


% ITERNAL FUCTION

function tftbcomp(path_dest)
%---------------------------------------------------------------
% compilation and displacement of the result to the directory
% 'directory'
%---------------------------------------------------------------

progs={'Cambifunb'
       'Ctfrdist'
       'Ctfrker'
       'Ctfrreas'
       'Ctfrrsp'
       'Ctfrstft'
       'Caf2tfr'
       'Cwindow'
       'Ctfrbj'
       'Ctfrbud'
       'Ctfrpwv'
       'Ctfrspwv'
       'Ctfrwv'
       'Ctfrzam'
       'Ctfrsp'
       'Ctfrgrd'
       'Ctfrridb'
       'Ctfrridh'
       'Ctfrridt'
       'Ctfrridbn'
       'Ctfrcw'
       'Ctfrri'
       'Ctfrpage'
       'Ctfrppage'
       'Ctfrmh'
       'Ctfrmhs'
       'Ctfrpmh'
       'Ctfrmmce'
       'Chtl'
       };
% number of programs to compile
nbre=size(progs);
nb_progs = nbre(1);
h=waitbar(0,'Compiling');
for i=1:nb_progs,
  waitbar(i/nb_progs,h);
  
  % compiles the C file
  name=char(progs(i)); 
  name=strcat(name,'.c');
  fprintf('%d ---> %s \n',i,name);
  name=fullfile('src',name);
  command='mex';
  command(end+1)=' ';
  command(end+1:end+length(name))=name;
  eval(command);

  % moves the mex file to its final destination
  name=char(progs(i)); 
  name=strcat(name,'.',mexext);
  copyfile(name,fullfile(path_dest,name));
  delete(name);
  
  % copies the help .m file in its final destination
  name=char(progs(i)); 
  name=strcat(name,'.m');
  copyfile(fullfile('hlpfiles',name),fullfile(path_dest,name));
end
copyfile(fullfile('hlpfiles','Contents.m'),fullfile(path_dest,'Contents.m'));
copyfile(fullfile('hlpfiles','Ctftbdemo.m'),fullfile(path_dest,'Ctftbdemo.m'));
close(h);
