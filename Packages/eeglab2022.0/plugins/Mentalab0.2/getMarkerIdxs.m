function [eeg_marker, orn_marker] = getMarkerIdxs(exg_timestamp, orn_timestamp, marker)
% GETMARKERS Gets the index of the ExG and ORN timestamps that correspond
% closest to the actual marker times
%   
% To do this, we find the two indices of the timestamp data that sit 
% either side of the marker time. That is, if in the exg_timestamp we have
% 1.1,1.2,1.3,1.4 and our marker timestamp is 1.25, we will get index 2 and
% 3. Then we see which timestamp is closest and return this value.

    eeg_marker = cell(height(marker), 2);
    orn_marker = cell(height(marker), 2);
    for i = 1:size(marker, 1)
        idx_eeg_above_marker = find(exg_timestamp > marker{i, 1}, 1);
        idx_eeg_below_marker = idx_eeg_above_marker - 1;
        
        diff_eeg_above = abs(exg_timestamp(idx_eeg_above_marker) - marker{i, 1});
        diff_eeg_below = abs(exg_timestamp(idx_eeg_below_marker) - marker{i, 1});

        %latency
        if (diff_eeg_above > diff_eeg_below)
            eeg_marker{i, 1} = idx_eeg_below_marker;
        else
            eeg_marker{i, 1} = idx_eeg_above_marker;
        end

        idx_orn_above_marker = find(orn_timestamp > marker{i, 1}, 1);
        idx_orn_below_marker = idx_orn_above_marker - 1;

        diff_orn_above = abs(orn_timestamp(idx_orn_above_marker) - marker{i, 1});
        diff_orn_below = abs(orn_timestamp(idx_orn_below_marker) - marker{i, 1});

        if (diff_orn_above > diff_orn_below)
            orn_marker{i, 1} = idx_orn_below_marker;
        else
            orn_marker{i, 1} = idx_orn_above_marker;
        end

        %label
        if (~isnumeric(marker{i, 2}))
            eeg_marker{i, 2} = char(marker{i, 2});
            orn_marker{i, 2} = char(marker{i, 2});
        else
            eeg_marker{i, 2} = marker{i, 2};
            orn_marker{i, 2} = marker{i, 2};
        end
    end
end
