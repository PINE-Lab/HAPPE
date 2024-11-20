    function [relpath pathstring] = find_relative_path(path_1,path_2)
        % Outputs the relative path from path_2 to path_1
        % So that path_1 = [path_2 filesep relative_path]

        % If they are the same, there is no change.
        if strcmp(path_1,path_2)
            relpath=cell(0);
            pathstring=path_1;
            return;
        end
        
        % The paths must be on the same drive/mount path
        x=filesep;
        if x=='\'
            x='\\';
        end
        
        Current=textscan(path_1,'%s','Delimiter',x);
        Current=Current{1};
        Path1=Current;
        Current=textscan(path_2,'%s','Delimiter',x);
        Current=Current{1};
        Path2=Current;
        for u=1:length(Path2)
            if u>length(Path1) || ~strcmp(Path1{u},Path2{u})
                break;
            end
        end
        paths_diverge_at=u;
        if paths_diverge_at==1
            error('Paths must contain at least one common directory.');
        end
        relpath=cell(0);
        for u=length(Path2):-1:paths_diverge_at
            relpath{end+1}='..';
        end
        for u=paths_diverge_at:length(Path1)
            relpath{end+1}=Path1{u};
        end
        
        pathstring='';
        start_from=1;
        tmpath=cd;
        for u=1:length(relpath)
            pathstring=[pathstring filesep relpath{u}];
            cd([path_2 filesep pathstring]);
            if strcmp(cd,path_2)
                start_from=u+1;
            end
        end
        if start_from<=length(relpath)
            relpath={relpath{start_from:end}};
        else
            relpath=cell(0);
        end
        cd(tmpath);
    end