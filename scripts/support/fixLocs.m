function locs = fixLocs(origLocs, happeDir)
    locs_1020 = readlocs([happeDir filesep 'Packages' filesep 'eeglab2024.0' ...
        filesep 'sample_locs' filesep 'Standard-10-20-Cap81.locs']) ;
    defaultChanIDs = {locs_1020.labels} ;
    locs = struct('labels', [], 'type', [], 'theta', [], 'radius', [], ...
        'X', [], 'Y', [], 'Z', [], 'sph_theta', [], 'sph_phi', [], ...
        'urchan', [], 'ref', []) ;
        % 'sph_radius', [], 'urchan', [], 'ref', []) ;

    for currChan = 1:size(origLocs, 2)
        try
            chanIndx = find(strcmp(defaultChanIDs, ...
                origLocs{currChan})) ;
            locs(currChan).labels = origLocs{currChan} ;
            locs(currChan).theta = locs_1020(chanIndx).theta ;
            locs(currChan).radius = locs_1020(chanIndx).radius ;
            locs(currChan).X = locs_1020(chanIndx).X ;
            locs(currChan).Y = locs_1020(chanIndx).Y ;
            locs(currChan).Z = locs_1020(chanIndx).Z ;
            locs(currChan).sph_theta = locs_1020(chanIndx).sph_theta ;
            locs(currChan).sph_phi = locs_1020(chanIndx).sph_phi ;
            % locs(currChan).sph_radius = locs_1020(chanIndx).sph_radius ;
            locs(currChan).type = locs_1020(chanIndx).type ; 
            locs(currChan).urchan = currChan ;
        catch ME
        end
    end
end