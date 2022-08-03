function out = byte2int24(input)
% This function gets 495 bytes (one data packet of EEG4) and converts it to
% voltage values

data = uint32(reshape(input, 3, size(input, 1) / 3)');
temp2 = int32(data(:,1) + data(:,2) * 2^8 + data(:,3) * 2^16);
out = temp2 - int32((temp2 >= 2^23) * (2^24));

end