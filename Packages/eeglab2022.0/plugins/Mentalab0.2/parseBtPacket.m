function [output, shownErr] = parseBtPacket(fid, shownErr)
% Parse the incoming data package
%   This function reads and parses one data package from the input stream.
%   
%   Input:
%       - fid: file/stream identifier
%   Output:
%       - output: structure of parsed data which depends on its type
%               Common fields:
%                   - cnt: counter of package
%                   - timestamp: time stamp of package
%               Type: Orientation
%                   - output.type:'end'
%                   - output.orn: Vector of Orientation data
%               Type: Environment
%                   - output.type: 'env'
%                   - output.temperature: in Celsius
%                   - output.light: in Lux
%                   - output.battery: in Volts
%               Type: Timestamp
%                   - output.type: 'ts'
%                   - output.timestamp: time stamp in seconds
%               Type: EEG4 or EEG8
%                   - output.type:'eeg4' or 'eeg8'
%                   - output.data: matrix of data (dimensions depend on the
%                   device type)
%               Type: end
%                   - output.type:'end'
%                   You get this output if the file/stream has ended or
%                   interrupted
%
%   Github page: https://github.com/Mentalab-hub/explorematlab/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = [];

EXG_UNIT = 1e-6;
TIMESTAMP_SCALE = 10000;
interruptWarning = 'Unexpected end of stream.';
endOfFileWarning = 'End of file.';
fletcherMismatchWarning = 'Fletcher mismatch!';
pidUnexpectedWarning = 'Unexpected package ID: ';

[pid, n] = fread(fid, 1, 'uint8');    % Read the package ID
if n == 0
    warning(endOfFileWarning);
    output.type = 'end';
    return; % No data in the stream/file
end

output.cnt = fread(fid, 1, 'uint8');                           % Counter of the package
payload = fread(fid, 1, 'uint16');                             % Number of bytes in the package
output.timestamp = fread(fid, 1, 'uint32') / TIMESTAMP_SCALE;  % Timestamp of the package in second

switch pid
    case 13                                                    % Orientation package
        output.type = 'orn';
        output.orn = fread(fid, (payload-8)/2, 'int16');
        if (numel(output.orn) < 9)
            uiwait(msgbox("The file you have requested is incomplete. I will attempt to fix missing data.","Error","error"));
            output.orn = [output.orn; zeros(9 - numel(output.orn), 1)];
            % I don't know if this is the right thing to do.
            % Essentially this seems to occur in the final item of a binary
            % file; it has been incompletely written. So I put
            % zeroes at the end to fill in the missing values...
        end
        output.orn = output.orn .* [0.061, 0.061, 0.061, 8.750, 8.750, 8.750, 1.52, 1.52, 1.52]';
    case {144, 146, 30, 62, 208, 210}                          % EEG package
        [temp, n] = fread(fid, (payload-8), 'uint8');
        if n < (payload-8) % check if the package terminates in between
            warning(interruptWarning);
            output.type = 'end';
            return;
        end
        if (pid == 144) || (pid == 208) % Specify the number of channel and reference voltage
            output.type = 'eeg4';
            nChan = 5; % 4 channels + 1 status
            vref = 2.4;
            nPacket = 33;
            temp = byte2int24(temp);
            [temp, shownErr] = reshapeWithWarning(temp, nChan, nPacket, shownErr);
            output.data = double(temp(2:end, :)) * vref / (2^23 - 1) / 6; % Calculate the real voltage value
        elseif (pid == 146) || (pid == 210)
            output.type = 'eeg8';
            nChan = 9; % 8 channels + 1 status
            vref = 2.4;
            nPacket = 16;
            temp = byte2int24(temp);
            [temp, shownErr] = reshapeWithWarning(temp, nChan, nPacket, shownErr);
            output.data = double(temp(2:end, :)) * vref / ( 2^23 - 1 ) / 6; % Calculate the real voltage value
        elseif pid == 30
            output.type = 'eeg8';
            nChan = 9; % 8 channels + 1 status
            vref = 4.5;
            nPacket = 16;
            temp = byte2int24(temp);
            [temp, shownErr] = reshapeWithWarning(temp, nChan, nPacket, shownErr);
            output.data = double(temp(2:end, :)) * vref / (2^23 - 1) / 6; % Calculate the real voltage value
        elseif pid == 62
            output.type = 'eeg8';
            nChan = 8;
            vref = 4.5;
            nPacket = 18;
            temp = byte2int24(temp);
            [temp, shownErr] = reshapeWithWarning(temp, nChan, nPacket, shownErr);
            output.data = double(temp) * vref / (2^23 - 1) / 6; % Calculate the real voltage value
        end
        output.data = round(output.data / EXG_UNIT, 2);
    case 194
        output.type = 'marker_event';
        output.code = fread(fid, 1, 'uint16');
    case 99
        output.type = 'dev_info';
        fw_str = num2str(fread(fid, 1, 'uint16'));
        output.fw_version = [fw_str(1) '.' fw_str(2) '.' fw_str(3)];
        output.data_rate = 16000 / (2 ^ fread(fid, 1, 'uint8'));
        output.adc_mask = dec2bin(fread(fid, 1, 'uint8'), 8);
    case {27, 19, 111, 192, 193, 195} % Not implemented / do nothing
        fread(fid, payload-8, 'uint8');
        output.type = 'unimplemented';
    otherwise
        warning([pidUnexpectedWarning pid])
        temp = fread(fid, payload-8, 'uint8'); % Read the payload
        output.type = 'end';
end


% Check the consistency of the Fletcher
[fletcher, n] = fread(fid, 4, 'uint8');
if n < 4
    warning(endOfFileWarning);
elseif ((pid ~= 27) && (fletcher(4) ~= 222))...
       || ((pid == 27) && (fletcher(4) ~= 255))
    disp(fletcher);
    warning(fletcherMismatchWarning)
    output.type = 'end';
end

end

function [temp, shownError] = reshapeWithWarning(temp, nChan, nPacket, shownError)
    if (abs(numel(temp) - nChan * nPacket) > 1) % #channels is wrong
        if (~shownError)
            uiwait(msgbox("The file you have requested is corrupt. I will attempt to fix it, however this may cause noticeable problems downstream.","Error","error"));
            shownError = true;
        end
        % channels are missing filled with zero
        temp = [temp; zeros(abs(nChan * nPacket - numel(temp)), 1)];
    end
    temp = reshape(temp, [nChan, nPacket]);
end
