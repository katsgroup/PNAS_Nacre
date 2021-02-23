function sampleRpart = ratioPart_IQ(fidr, spatial, spectral, framesUsed, whiteAvg, darkAvg, dark4whiteAvg)

V = fread(fidr,spatial*spectral*framesUsed,'uint16=>uint16');
sample = single(reshape(V,spatial,spectral,framesUsed));
sample = permute(sample,[3,1,2]);
whiteFull = single(repmat(whiteAvg,framesUsed,1));
darkFull = single(repmat(darkAvg,framesUsed,1));
dark4whiteFull = single(repmat(dark4whiteAvg,framesUsed,1));

if size(sample)~=size(whiteFull)
    uiwait(msgbox('The sample and references are not the same number of spatial points. Cannot ratio sample.','Error','modal'));
    return;
end
sampleRpart = (sample - darkFull)./(whiteFull - dark4whiteFull);