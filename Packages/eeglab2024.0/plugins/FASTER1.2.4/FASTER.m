function FASTER(option_wrapper)

% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

startDir=option_wrapper.options.file_options.folder_name;
outDir=option_wrapper.options.file_options.output_folder_name;
chan_locs=option_wrapper.options.file_options.channel_locations;
is_bdf=option_wrapper.options.file_options.is_bdf==1;
using_ALLEEG=option_wrapper.options.file_options.using_ALLEEG;
resume=option_wrapper.options.file_options.resume;

if ~using_ALLEEG
    %[jpathname jf1 jf2 jf3]=fileparts(option_wrapper.options.job_filename);
    [jpathname jf1 jf2]=fileparts(option_wrapper.options.job_filename);
    Qname=[jf1 jf2 '_ProcQ.eegQ'];
else
    Qname='ProcQ.eegQ';
end

is_processing=exist([startDir filesep Qname],'file');

if ((~ischar(chan_locs) || ~exist(chan_locs,'file')) && is_bdf) && ~is_processing
    fprintf('Invalid channel location file.\n');
    return;
end

if (~using_ALLEEG)
    can_distribute=1;
    if ~is_processing
        if (resume)
            if ~isempty(option_wrapper.options.file_options.plist)
                plist=option_wrapper.options.file_options.plist;
                nlist=option_wrapper.options.file_options.nlist;
                oplist=option_wrapper.options.file_options.oplist;
            else
                if (is_bdf)
                    [plist nlist] = extsearchc(startDir,'.bdf',0);
                else
                    % Assume .set if not .bdf
                    % Later versions may support other file formats
                    [plist nlist] = extsearchc(startDir,'.set',0);
                end
                oplist=cell(size(plist));

                x=true(size(plist));
                for i=1:length(plist)
                    if length(plist{i})>12
                        if strcmp(plist{i}(end-11:end),'Intermediate')
                            x(i)=0;
                        end
                    end
                    if isempty(outDir)
                        oplist{i}=plist{i};
                    else
                        oplist{i}=outDir;
                        dir_structure=get_dir_structure(plist{i},startDir);
                        for k=1:length(dir_structure)
                            oplist{i}=[filepath filesep dir_structure{k}];
                            if ~exist(oplist{i},'dir')
                                mkdir(oplist{i});
                            end
                        end
                    end
                end
                plist={plist{x}};
                nlist={nlist{x}};
                oplist={oplist{x}};

                option_wrapper.options.file_options.plist=plist;
                option_wrapper.options.file_options.nlist=nlist;
                option_wrapper.options.file_options.oplist=oplist;
            end
        else
            if (is_bdf)
                [plist nlist] = extsearchc(startDir,'.bdf',0);
            else
                % Assume .set if not .bdf
                % Later versions may support other file formats
                [plist nlist] = extsearchc(startDir,'.set',0);
            end
            oplist=cell(size(plist));
            x=true(size(plist));
            
            filepath=option_wrapper.options.file_options.output_folder_name;
            
            for i=1:length(plist)
                if length(plist{i})>12
                    if strcmp(plist{i}(end-11:end),'Intermediate')
                        x(i)=0;
                    end
                end
                if isempty(outDir)
                    oplist{i}=plist{i};
                else
                    oplist{i}=outDir;
                    dir_structure=get_dir_structure(plist{i},startDir);
                    for k=1:length(dir_structure)
                        oplist{i}=[filepath filesep dir_structure{k}];
                        if ~exist(oplist{i},'dir')
                            mkdir(oplist{i});
                        end
                    end
                end
            end
            plist={plist{x}};
            nlist={nlist{x}};
            oplist={oplist{x}};
            if (option_wrapper.options.file_options.make_subdirectories)
                x = findrepeats(plist);
                for i = 1:length(x)
                    if isempty(outDir)
                        try
                            mkdir(plist{x(i)},nlist{x(i)}(1:end-4));
                            movefile([plist{x(i)} filesep nlist{x(i)}],[plist{x(i)} filesep nlist{x(i)}(1:end-4)],'f');
                            if exist([plist{x(i)} filesep nlist{x(i)}(1:end-4) '.fdt'],'file') && ~is_bdf
                                movefile([plist{x(i)} filesep nlist{x(i)}(1:end-4) '.fdt'],[plist{x(i)} filesep nlist{x(i)}(1:end-4)],'f');
                            end
                            plist{x(i)} = [plist{x(i)} filesep nlist{x(i)}(1:end-4)];
                            oplist{x(i)} = plist{x(i)};
                        catch
                            error('Error in organising files in %s\n',plist{x(i)});
                            return;
                        end
                    else
                        try
                            if ~exist([outDir filesep nlist{x(i)}(1:end-4)],'dir')
                                mkdir(outDir,nlist{x(i)}(1:end-4));
                            end
                            oplist{x(i)} = [outDir filesep nlist{x(i)}(1:end-4)];
                        catch
                            error('Error in organising files in %s\n',plist{x(i)});
                            return;
                        end
                    end
                end
            end
            option_wrapper.options.file_options.plist=plist;
            option_wrapper.options.file_options.nlist=nlist;
            option_wrapper.options.file_options.oplist=oplist;
        end

        % Changed the below from copyfile to direct write due to some weird
        % permission issue that was probably only on one test computer.
        % Works now anyway!
        if ~isempty(chan_locs)
            [clpathstr, clname, clext] = fileparts(chan_locs);
            fid=0;
            if ~isempty(outDir) && ~strcmp(outDir,clpathstr)
                %copyfile(chan_locs,outDir,'f');
                fid=fopen([outDir filesep clname clext],'w');
            elseif ~strcmp(startDir,clpathstr)
                %copyfile(chan_locs,startDir,'f');
                fid=fopen([startDir filesep clname clext],'w');
            end
            if fid>0
            fid2=fopen(chan_locs,'r');
            while ~feof(fid2)
                fprintf(fid,'%s',fgets(fid2));
            end
            fclose(fid);
            fclose(fid2);
            end
        end

        % Add stuff to save in the [startDir filesep 'ProcQ.eegQ'] file
        Q.plist=plist;
        Q.nlist=nlist;
        Q.plist_rel=cell(0);
        Q.oplist_rel=cell(0);
        for v=1:length(plist)
            Q.plist_rel{v}=find_relative_path(plist{v},startDir);
            Q.oplist_rel{v}=find_relative_path(oplist{v},startDir);
        end
        Q.outDir_rel=cell(0);
        if ~isempty(outDir)
            Q.outDir_rel=find_relative_path(outDir,startDir);
        end
        my_comp_num=1;
        Q.comp_nums=1;
        Q.finished=0;
        Q.processed=zeros(size(Q.plist));
        Q.errors=zeros(size(Q.plist));
        if resume
            Q.next_file=option_wrapper.options.file_options.current_file_num;
        else
            Q.next_file=1;
        end
        save([startDir filesep Qname],'Q');
    else
        my_queue_file=[];
        wait_and_lock();
        L=load([startDir filesep Qname],'-mat');
        Q=L.Q;
        clear L;
        % Recreate the path list from the relative path list, which are
        plist=cell(size(Q.plist_rel));
        oplist=plist;
        for v=1:length(Q.plist_rel)
            plist{v} = make_relative_path(Q.plist_rel{v},startDir);
            oplist{v} = make_relative_path(Q.oplist_rel{v},startDir);
        end
        %plist=Q.plist;
        nlist=Q.nlist;
        if ~isempty(Q.outDir_rel)
            outDir = make_relative_path(Q.outDir_rel,startDir);
        end
        %fprintf('\n\nOutput directory is: %s\n\n',outDir);
        my_comp_num=max(Q.comp_nums)+1;
        Q.comp_nums=[Q.comp_nums my_comp_num];
        save_and_unlock();
        [clpathstr, clname, clext] = fileparts(chan_locs);
        if ~isempty(outDir)
            option_wrapper.options.file_options.channel_locations=[outDir filesep clname clext];
        else
            option_wrapper.options.file_options.channel_locations=[startDir filesep clname clext];
        end
    end
else
    can_distribute=0;
    if is_processing==1
        error('Can''t join a running jobfile if opened from EEGLAB. If this file is not being processed, reload it and reset the queue.');

        % The below won't run, I may reinstate it at some point, but it's
        % difficult to co-ordinate and may be pointless
        c=clock;
        months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
        new_jobname=[sprintf('%d/%s/%d__%d-%d-%d___',c(3),months{c(2)},c(1),c(4),c(5),round(c(6))) option_wrapper.options.job_filename];
        warning(sprintf('This jobfile is running on another machine.\nA new jobfile will be created: %s',new_jobname)); %#ok<WNTAG,SPWRN>
        option_wrapper.option.job_filename=new_jobname;
        return;
    end
    startDir=cd;
    plist=zeros(evalin('base','size(ALLEEG);'));
    %nlist=plist;
    for v=1:length(plist)
        plist(v)=evalin('base',sprintf('isempty(ALLEEG(%d).data);',v));
    end
    ALLEEG_to_do=find(plist~=1);
    plist=cell(length(ALLEEG_to_do));
    for v=1:length(plist)
        plist{v}=ALLEEG_to_do(v);
        oplist=cell(size(plist));
        if ~isempty(outDir)
            oplist{v}=outDir;
        else
            oplist{v}=evalin('base',sprintf('(ALLEEG(%d).filepath);',plist{v}));
            if isempty(oplist{v})
                oplist{v}=cd;
            end
        end
    end
    nlist=plist;
    Q.plist=plist;
    Q.nlist=nlist;
    Q.plist_rel=plist;
    my_comp_num=1;
    Q.comp_nums=1;
    Q.finished=0;
    Q.processed=zeros(size(Q.plist));
    Q.errors=zeros(size(Q.plist));
    %Q.next_file=ALLEEG_to_do(1);
    Q.next_file=1;
    save([startDir filesep Qname],'Q');
    option_wrapper.options.file_options.plist=plist;
    option_wrapper.options.file_options.nlist=nlist;
    option_wrapper.options.file_options.oplist=oplist;
end

if (~exist([startDir filesep 'Processing'],'dir'))
    mkdir([startDir filesep 'Processing']);
end
if (~exist([startDir filesep 'Queue'],'dir'))
    mkdir([startDir filesep 'Queue']);
end

all_errors=cell(0);

doing_distrib=0;
is_processing=1;
error_indices=zeros(size(plist));
first_file=1;
had_error=0;
my_proc_file=[];
make_processing_file();
EEG_state=[];
while (1)
    %%%%%%%%%%%%%%%%%%%%%
    % Before processing %
    %%%%%%%%%%%%%%%%%%%%%

    % The queue system is there to ensure that multiple computers are not
    % reading the options file simultaneously and so trying to process the
    % same file! Also to make sure all computers have finished processing
    % before the grand average is made.

    % Puts the current computer number in the queue
    % and waits until other computers are done updating
    % the queue file. Then opens the queue file, reads the
    % current state to find the next file to process, then
    % updates that for the next computer to read. If there
    % is no next file, it marks the processing as finished.
    my_queue_file=[];
    wait_and_lock();
    L=load([startDir filesep Qname],'-mat');
    Q=L.Q;
    if ~first_file
        if ~had_error
            Q.processed(current_file)=1;
        else
            Q.errors(current_file)=1;
            if option_wrapper.debug && exist('m','var')
                if isempty(outDir)
                    if exist([startDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'file')
                        L=load([startDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'-mat');
                        all_errors=L.all_errors;
                    end
                    all_errors{end+1,1}=m;
                    all_errors{end,2}=EEG_state;
                    save([startDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'all_errors','-mat');
                    error_indices(current_file)=size(all_errors,1);
                else
                    if exist([outDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'file')
                        L=load([outDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'-mat');
                        all_errors=L.all_errors;
                    end
                    all_errors{end+1,1}=m;
                    all_errors{end,2}=EEG_state;
                    save([outDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'all_errors','-mat');
                    error_indices(current_file)=size(all_errors,1);
                end
            end
        end
    end
    had_error=0;
    %if (~using_ALLEEG && Q.next_file>length(Q.plist)) || (using_ALLEEG && Q.next_file>max(ALLEEG_to_do))
    if Q.next_file>length(Q.plist)
        Q.finished=1;
    end
    if Q.finished
        Q.comp_nums=setdiff(Q.comp_nums,my_comp_num);
        delete_processing_file();
        save_and_unlock();
        break;
    end

    current_file=Q.next_file;
    Q.next_file=Q.next_file+1;

    option_wrapper.options.file_options.current_file = [plist{current_file} filesep nlist{current_file}];
    option_wrapper.options.file_options.current_file_num=current_file;
    if (~isempty(option_wrapper.options.job_filename))
        save(option_wrapper.options.job_filename,'option_wrapper','-mat');
    end

    save_and_unlock();

    %%%%%%%%%%%%%%
    % Processing %
    %%%%%%%%%%%%%%
    if (~isempty(option_wrapper.options.file_options.searchstring))
        searchstring2=option_wrapper.options.file_options.searchstring;
    else
        searchstring2=nlist{current_file};
    end
    if ((using_ALLEEG && ~evalin('base',sprintf('isempty(ALLEEG(%d).data);',plist{current_file}))) || ~isempty(strfind(nlist{current_file},searchstring2)) )
        tic
        fprintf('******************\n');
        fprintf('* File %.3d / %.3d *\n',current_file,length(nlist));
        fprintf('******************\n');

        if (~using_ALLEEG)
            %log_file = fopen([plist{current_file} filesep option_wrapper.options.file_options.file_prefix nlist{current_file}(1:end-4) '.log'],'a+');
            log_file = fopen([oplist{current_file} filesep option_wrapper.options.file_options.file_prefix nlist{current_file}(1:end-4) '.log'],'a+');
        else
            log_file = fopen([oplist{current_file} filesep option_wrapper.options.file_options.file_prefix sprintf('FASTER_ALLEEG(%d).log',ALLEEG_to_do(current_file))],'a+');
        end

        try
            FASTER_process(option_wrapper,log_file);
        catch
            m=lasterror;
            fprintf('\nError - %s.\n',m.message);
            try fclose(log_file); catch; end;
            had_error=1;
        end
    else
        fprintf('Skipped file.\n');
    end

    % After processing
    first_file=0;
end

%%%%%%%%%%%%%%%%%%%%%%
% Queueing Functions %
%%%%%%%%%%%%%%%%%%%%%%
    function wait_and_lock()
        % Find the last queue file, make one for this computer, then wait
        % until the previous one is deleted
        D=dir([startDir filesep 'Queue']);
        N={D(:).name};
        N=setdiff(N,{'.','..'});
        N=str2double(N);
        if isempty(N)
            N=0;
        end
        next_queue_num=max(N)+1;
        my_queue_file=[startDir filesep 'Queue' filesep sprintf('%d',next_queue_num)];
        fid=fopen(my_queue_file,'w');
        fclose(fid);
        prev_queue_file=[startDir filesep 'Queue' filesep sprintf('%d',max(N))];
        said=0;
        while max(N)>0 && exist(prev_queue_file,'file')
            if ~said
                fprintf('Waiting for removal of Queue file %s\n',prev_queue_file);
                said=1;
            end
            pause(1);
        end
        assignin('caller','my_queue_file',my_queue_file);
    end
    function save_and_unlock()
        save([startDir filesep Qname],'Q');
        pause(1);
        delete(my_queue_file);
        assignin('caller','my_queue_file',[]);
    end
    function make_processing_file()
        my_proc_file=[startDir filesep 'Processing' filesep sprintf('%d',my_comp_num)];
        fid=fopen(my_proc_file,'w');
        fclose(fid);
        assignin('caller','my_proc_file',my_proc_file);
    end
    function delete_processing_file()
        delete(my_proc_file);
        assignin('caller','my_proc_file',[]);
    end

%%%%%%%%%%%%%%%%%%%
% Post processing %
%%%%%%%%%%%%%%%%%%%

if using_ALLEEG
    evalin('base','EEG=ALLEEG(CURRENTSET);');
    evalin('base','eeglab redraw;');
end

D=dir([startDir filesep 'Processing']);
if length(D)>2
    fprintf('*******************\n');
    fprintf('* FASTER Finished *\n');
    fprintf('*******************\n');
    fprintf('Finished processing all my files. The last computer to finish processing will make the grand average.\n');
    return;
end

if isempty(outDir)
    top_log = fopen([startDir filesep option_wrapper.options.file_options.file_prefix 'FASTER.log'],'a');
    if exist([startDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'file')
        L=load([startDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'-mat');
        all_errors=L.all_errors;
    end
else
    top_log = fopen([outDir filesep option_wrapper.options.file_options.file_prefix 'FASTER.log'],'a');
    if exist([outDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'file')
        L=load([outDir filesep option_wrapper.options.file_options.file_prefix 'FASTER_errors.mat'],'-mat');
        all_errors=L.all_errors;
    end
end

c=clock;
months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
fprintf(top_log,'\n%d/%s/%d %d:%d:%d\n',c(3),months{c(2)},c(1),c(4),c(5),round(c(6)));

for v=1:length(plist)
    if ~using_ALLEEG
        fprintf(top_log,'%s%s%s:\n',plist{v},filesep,nlist{v});
    else
        fprintf(top_log,'ALLEEG(%d):\n',v);
        if Q.processed(v)
            fprintf(top_log,'Processed successfully.\n');
        elseif Q.errors(v)
            fprintf(top_log,'Error: %s. Load all_errors.mat and investigate all_errors{%d,1} for more info.\n',all_errors{error_indices(v),1}.message,error_indices(v));
        else
            fprintf(top_log,'Skipped due to filename filter or empty ALLEEG dataset.\n');
        end
        fprintf(top_log,'\n');
    end
end

delete([startDir filesep Qname]);
if length(dir([startDir filesep 'Processing']))==2
    rmdir([startDir filesep 'Processing']);
end
if length(dir([startDir filesep 'Queue']))==2
    rmdir([startDir filesep 'Queue']);
end

%%%%%%%%%%%%%%%%%
% Grand Average %
%%%%%%%%%%%%%%%%%

FASTER_grandaverage(startDir,option_wrapper,all_errors,top_log,plist,nlist);

%%%%%%%%%%%%%
% All done! %
%%%%%%%%%%%%%

fprintf('*******************\n');
fprintf('* FASTER Finished *\n');
fprintf('*  %.3d processed  *\n',sum(Q.processed));
fprintf('*   %.3d errors    *\n',sum(Q.errors));
fprintf('*   %.3d skipped   *\n',length(plist)-sum(Q.processed)-sum(Q.errors));
fprintf('*******************\n');
fprintf(top_log,'\nFinished. %d processed, %d errors, %d skipped.\n',sum(Q.processed),sum(Q.errors),length(plist)-sum(Q.processed)-sum(Q.errors));
fclose(top_log);

    function out_chan_locs = check_chan_locs(chan_locs,num_chans,num_exts)
        fid = fopen(chan_locs,'r');
        if fid==-1
            out_chan_locs = -1;
            return;
        end
        count=0;
        x='';
        while ~feof(fid)
            x = fgetl(fid);
            if ~isempty(x)
                count = count+1;
            end
        end
        if count < num_chans + num_exts
            out_chan_locs = -1;
            fclose(fid);
            return;
        end
        if count == num_chans + num_exts
            out_chan_locs = chan_locs;
            fclose(fid);
            return;
        end
        frewind(fid);
        fid2 = fopen(['tempchanloc' chan_locs(end-4:end)],'w');
        count2 = 0;
        while ~feof(fid)
            x = fgets(fid);
            if (count2 < num_chans || count - count2 <= num_exts)
                fprintf(fid2,'%s',x);
            end
            count2 = count2 + 1;
        end
        fclose(fid);
        fclose(fid2);
        out_chan_locs = ['tempchanloc' chan_locs(end-4:end)];
    end

    function indices = findrepeats(input)
        indices=zeros(size(input));
        for u=1:length(input)
            indices(u)=sum(strcmp(input{u},input));
        end
        indices=find(indices>1);
    end

end