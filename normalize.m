function result = normalize(spectrum, maxIntensity)
sz = size(spectrum);
result = zeros(sz);
max = -1;
for i = 1:sz(2)
    if spectrum(i) > max
        max = spectrum(i);
    end
end
factor = maxIntensity/max;
for i = 1:sz(2)
     result(i) = spectrum(i)*factor;
end