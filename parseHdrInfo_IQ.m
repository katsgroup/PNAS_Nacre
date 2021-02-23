%% parseHdrInfo
%
% [wavelengths, spatial, frames, spectral, tint, settings ] = parseHdrInfo(path,name)

%% Begin Code
function [ wavelengths, spatial, frames, spectral, tint, settings ] = parseHdrInfo_IQ( path, name )
                %Enter a path (e.g. 'M:\Engineering\...') and .hdr filename
f = fullfile(path, name);       
fp = fopen(f, 'r');             %Opens file as read-only

spatial =0;                     %Initialize spatial and spectral
spectral=0;
waveFlag = 0;                   %Flag changes to 1 when wavelengths vector is collected
wavelengths = zeros(1,spectral);
settings = zeros(1,spectral);

line = fgetl(fp);               %Sets variable line to the first line (and then the next, the next) of the file
while(~feof(fp)) && (waveFlag == 0)   %While not end of file and spatial or waveflag is 0
    ind = strfind(line, 'samples');
    if (ind > 0)
        spatial = str2double(line((strfind(line, '=')+1):size(line,2)));
    else
        ind=strfind(line,'lines');
        if (ind > 0)
            frames = str2double(line((strfind(line,'=')+1):size(line,2)));
        else
            ind = strfind(line, 'bands');
            if (ind == 1)
                spectral = str2double(line((strfind(line, '=')+1):size(line, 2)));
            else
                ind = strfind(line, 'tint');
                if (ind == 1)
                    tint = str2double(line((strfind(line, '=')+1):size(line, 2)));
                else
                    ind = strfind(line, 'wavelength');     %If line has the word Wavelength, ind = character # of location of 'Wavelength' (1, b/c line starts w/"Wavelength")
                    if(ind > 0 )
                        line=fgetl(fp);
                        for k = 1:spectral
                            wavelengths(1,k) = str2double(line(1:size(line,2)));
                            %Str2double converts argument text > vector
                            %Use characters from beginning to end of line (2nd dim.)
                            line=fgetl(fp);
                        end
%                         waveFlag = 1;
                    else
                        ind = strfind(line, 'Band Intensities');     %If line has the word Wavelength, ind = character # of location of 'Wavelength' (1, b/c line starts w/"Wavelength")
                        if(ind > 0 )
                            line=fgetl(fp);
                            for k = 1:spectral
                                settings(1,k) = str2double(line(1:size(line,2)));
                                %Str2double converts argument text > vector
                                %Use characters from beginning to end of line (2nd dim.)
                                line=fgetl(fp);
                            end
                            waveFlag = 1;
                        end
                    end
                end
            end
        end
    end
    line=fgetl(fp);
end
fclose(fp);

wavelengths = sort(wavelengths);
    
                       

                    
    
    
    
    