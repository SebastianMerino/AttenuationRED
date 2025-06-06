function bMode = getBmode(rf,fs,bw)
N = size(rf,1);
rfExt = wextend("addrow","sym",rf,N,"d");
[b,a] = butter(2,bw/fs*2,"bandpass");
refFilt = filtfilt(b,a,rfExt);
bMode = db(hilbert(refFilt));
bMode = bMode(1:N,:);
bMode = bMode - max(bMode(:));


end