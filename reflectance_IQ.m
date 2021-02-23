function R = reflectance_IQ(varargin)
%REFLECTANCE
% R = reflectance(whiteAvg,darkAvg,dark4whiteAvg,path,filename);
% Optional 6th argument: range of frames to use e.g. fr1:fr2

whiteAvg = varargin{1};
darkAvg = varargin{2};
dark4whiteAvg = varargin{3};
path = varargin{4};
filename = varargin{5};
skip = 0;
[~, spatial, frames, spectral] = parseHdrInfo_IQ(path, strrep(filename,'.raw','.hdr'));
if nargin>5
    frames = varargin{6};
    fr1 = frames(1);
    fr2 = frames(end);
    frames = fr2-fr1+1;
    skip = 1;
end

if skip
    pointsToSkip = spatial*spectral*(fr1-1)*2;
else
    pointsToSkip = 0;
end
framesUsed = min(25,frames);
numParts = floor(frames/framesUsed);
extraFr = frames-numParts*framesUsed;

fidr = fopen(fullfile(path,filename),'r');
fseek(fidr, pointsToSkip, 'bof');
R = single(zeros(frames,spatial,spectral));
for i=1:numParts
    R((i-1)*framesUsed+1:i*framesUsed,:,:) = ratioPart(fidr, spatial,...
        spectral, framesUsed, whiteAvg, darkAvg, dark4whiteAvg);
end
if extraFr>0
    R(framesUsed*numParts+1:frames,:,:) = ratioPart(fidr, spatial,...
        spectral, extraFr, whiteAvg, darkAvg, dark4whiteAvg);
end
fclose(fidr);