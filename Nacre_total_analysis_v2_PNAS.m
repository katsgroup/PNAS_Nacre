clear all;
close all;

%This code allows for taking a raw hyperspectral image of nacre, selecting
%the region of interest, then performs the fitting of thicknesses and
%standard deviations. Afterwards a sweeping analysis is done to create a
%final statistical representation of the data. The user is prompted to
%select the appropriate files for: hyperspectral data, header, white
%reference, simulated spectra.
%
%This file combines several scripts into one.

%09-09-2020: code adjusted in scrubbing section to only use a subset of the
%confidence map defined by the user's region of interest to determine the
%scrubbing threshold. 
%By Jad Salman
%% Extraction of Reflectances from hyperspectral data

% This is based on the Nacre_Hyperspectral_specim_IQ.m code. All files
% should be in present working directory.
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',3);

[file,path_sample] = uigetfile('*.raw','Select sample .raw file'); %Select .raw file of hyperspectral image.
[header_file,path_hdr] = uigetfile('*.hdr','Select sample .hdr file'); %Select .hdr file of hyperspectral image
[whiteFile,path_white] = uigetfile('*.raw','Select white .raw file'); %Select white reference .raw file
[darkForWhiteFile,path_darkforwhite] = uigetfile('DARK*.*raw','Select white DARK .raw file'); %Select .raw dark reference for the WHITE file
[darkForSampleFile,path_darkforsample] = uigetfile('DARK*.*raw','Select sample DARK .raw file'); %Select .raw the dark reference for the SAMPLE file
[wavelengths, spatial, frames, spectral, tint, settings ] = parseHdrInfo_IQ(path_sample,header_file); %Extract wavelength [nm] and other params.
%% Read in image file
img = reshapeImage_IQ(path_sample,file);
white = reshapeImage_IQ(path_white,whiteFile);
dark = reshapeImage_IQ(path_darkforwhite,darkForWhiteFile);
dark_sample = reshapeImage_IQ(path_darkforsample,darkForSampleFile);

wavlen = 60; %Define wavelength index to plot raw image
%Plot at a single wavelength and select region to consider as white
%reference white reference
figure
imagesc(white(:,:,wavlen)); %View single wavelength of white reference
title(['Raw white ref image @ \lambda = ',num2str(wavelengths(wavlen)),' nm']);
xlabel('Pixel number');
ylabel('Frame number');
axis 'equal'

[x_w,y_w] = getpts; %Grab region in which to average white reference
x_w = round(x_w); %Round to nearest integer
y_w = round(y_w);

%% Convert to reflectance

whiteAvg = zeros(spatial,spectral); %initialize average row of white pixels
darkavg = whiteAvg; %initialize average row of dark pixels
darkavg_img = whiteAvg;
img_dark_sub = zeros(frames,spatial,spectral); %initialize image matrix - dark pixels
white_sub = zeros(size(white,1),size(white,2),size(white,3));
whiteavg_sub = whiteAvg;

% Loop to find the average dark and white references. Subtracts the average
% dark from the image data.
for i = 1:spatial
    for j = 1:spectral
        darkavg(i,j) = mean(dark(:,i,j));
        darkavg_img(i,j) = mean(dark_sample(:,i,j));
        whiteAvg(i,j) = mean(mean(white(y_w(1):y_w(2),x_w(1):x_w(2),j)));
        whiteavg_sub(i,j) = whiteAvg(i,j)-darkavg(i,j);
        img_dark_sub(:,i,j) = img(:,i,j)-darkavg_img(i,j);
    end
end
%clear the raw image files to free up memory
clear dark
clear white
clear img

%Optional: Plot reflectance map at a single wavelength defined by wavelen
figure
imagesc(img_dark_sub(:,:,wavlen))
axis 'equal'
title(['Sample-dark image @ \lambda = ',num2str(wavelengths(wavlen)),' nm']);
xlabel('Pixel number');
ylabel('Frame number');

%Reflectance calculation. Normalize value in image to white reference
%R is normalized for each frame to the same row of white reference pixels
R = zeros(size(img_dark_sub));
for k = 1:frames
    for i = 1:spatial
        for j = 1:spectral
            R(k,i,j) = img_dark_sub(k,i,j)./whiteavg_sub(i,j);
        end
    end
end
%plot reflection 2D image at single wavelength
figure
imagesc(R(:,:,wavlen))
axis 'equal'
title(['Norm. Refl @ \lambda = ',num2str(wavelengths(wavlen)),' nm']);
xlabel('Pixel number');
ylabel('Frame number');

%% Extracting spatial region of hyperspectral data cube
%Click two points to define an area to grab spectra for each pixel
[x,y] = getpts;
x = round(x);
y = round(y);
xy = zeros(abs(x(2)-x(1))*abs(y(1)-y(2)),2);
count = 0;
for i = 1:abs(x(2)-x(1))
    for j = 1:abs(y(2)-y(1))
        count = count+1;
        xy(count,:) = [x(1)+i-1,y(1)+j-1];
    end
end

% OPTIONAL: Plot individual spectra for each pixel in selected area
% close all
% figure
% test_sum = 0;
% for i = 1:length(xy)
%     test = R(xy(i,2),xy(i,1),:);
% %     peak = max(test);
% %     test_norm = test./peak; %normalize all reflectance maximums to 1
%     plot(wavelengths, smooth(test(:),1))
%     hold on
% end

%Save matrix of pixel position points used to extract reflection data
pixel_matx = reshape(xy(:,1),abs(y(1)-y(2)),abs(x(1)-x(2)));
pixel_maty = reshape(xy(:,2),abs(y(1)-y(2)),abs(x(1)-x(2)));

%Save only R spectra for selected pixel positions
R_mat = R(y(1):1:y(2),x(1):1:x(2),:);
lambda_meas = wavelengths;

%% Part B: fitting of nacre data
%Fit hyperspectral data and map thickness to points across the
%hyperspectral map.
%
%This uses the mean-square error (MSE) of the derivatives of the normalized
%measured reflectance spectra and the normalized calculated reflectance
%spectra. Then a subset of each best fit point is refit using the
%nearest-neighbor spectra by performing a MSE calculation of the direct
%reflectance spectra.

%% Load data and set thickness ranges
%A precalculated set of simulated spectra is needed to load here. This is
%an extensive set of spectra calculated from transfer-matrix method
%combined with Monte-carlo method to simulate the response of the desired
%nacre. See code titled "____" for how this set of data is calculated.

%Step size for simulation spectra set based on step size used in the
%transfer-matrix calculations. Here, we defined a range and step size of
%150-600 nm in tablet thickness with 5 nm increments. The standard
%deviation was 5-70 nm with 5 nm increments.
t_nacre = 150:5:600; %Mean thickness range
sig_nacre = 5:5:70; %Mean standard deviation range

t_dim = length(t_nacre); % dimension of thickness fitting vector
sig_dim =length(sig_nacre); % dimension of standard dev fitting vector

%Select the appropriate simulated spectra
cd('Simulated Nacre R_out files'); % Point to directory where simulation spectra are located
load(uigetfile('*','Select R_out file'))
cd ..
%% Condition data
lambda_test = wavelengths;
R_out_trunc =R_out((10:end),:,:); %R_out is the simulated reflectances. %OPTIONAL: Truncate the first 10 data points if noisey
%R_out_trunc = R_out_trunc./max(R_out_trunc); %Re-normalize the simulated data to truncated peak
confidence = zeros(size(R_mat,1),size(R_mat,2)); %Initialize confidence value (peak intensity)/mse
% confidence2 = zeros(size(R_mat,1),size(R_mat,2)); %Alternative confidence metric: (mean intensity)/mse
% confidence3 = zeros(size(R_mat,1),size(R_mat,2)); %Alternative confidence metric: 1/mse

%Set an intensity threshold for spectra and eliminate all other spectra.
%This is OPTIONAL if you wish to premptively scrub spectra ahead of fitting
threshold = 0; %Reflectance peak intensity
tic
for k = 1:size(R_mat,1)
    for i = 1:size(R_mat,2)
        R_mat(k,i,:) = smooth(R_mat(k,i,:),1); %No smoothing
        if max(R_mat(k,i,19:end))< threshold
            R_mat(k,i,:) = nan;
        end        
    end
end
toc
% R_mat_smooth = smoothdata(R_mat,3)

%Measured spectra conditioning
sample = 1; % set down sample number based on simulated data sample rate
% lambda_down = downsample(lambda_test,sample); %wavelength downsample
lambda_down = lambda_test; %wavelength
R_meas_down = zeros(size(R_mat,1),size(R_mat,2),length(lambda_down(10:end))); %%Truncate the first 10 data points. CHANGE IF NUMBER OF TRUNCATED VALUES IS DIFFERENT
R_meas_trunc = zeros(size(R_mat,1),size(R_mat,2),length(lambda_down(10:end))); %No down sampling here. %Truncate the first 10 data points. CHANGE IF NUMBER OF TRUNCATED VALUES IS DIFFERENT

%Downsample, normalize, and truncate all reflectance data with respect to
%wavelength
lambda_trunc = lambda_down(10:end); % %Truncate the first 10 data points. CHANGE IF NUMBER OF TRUNCATED VALUES IS DIFFERENT


%figure %OPTIONAL: uncomment to show plots
for i= 1:size(R_mat,1)
    for j = 1:size(R_mat,2)
%         R_meas_down(i,j,:) = downsample(R_mat(i,j,10:end),sample);
        R_meas_down(i,j,:) = R_mat(i,j,10:end);
        test = R_meas_down(i,j,1:end); % pick single reflectance spectrum
        confidence(i,j) = max(test);
       % confidence2(i,j) = mean(test); %Different confidence metric
        norm_test = test(:)./max(test); %Normalize to peak
        R_meas_trunc(i,j,:) = norm_test; % Save truncated and normalized measured R
%         plot(lambda_trunc,norm_test); % Optional: uncomment to show plots
%         ylim([0,1])
%         hold on
        
    end
end

% Smooth truncated measured and simulated data for doing derivative of spectra
R_meas_trunc_smooth = smoothdata(R_meas_trunc,'movmean',10);
R_out_trunc_smooth = smoothdata(R_out_trunc,'movmean',10);

%% Fit calculated to measure R and map thicknesses

%Define dimensions of row, column, and page based on hyperspectral images
row = size(R_meas_trunc_smooth,1);
column = size(R_meas_trunc_smooth,2);
page = size(R_meas_trunc_smooth,3);

%Create cell matrix to handel best fit index for plotting
best_index = cell(row*column,1);


%Reshape 3-dimensional data sets to be 2 dimensions with the columns
%along the wavelength. Takes a m x n x p matrix and makes it (m*n) x p.
%Matrix values with _r are reshaped
R_meas_trunc_smooth_r = reshape(R_meas_trunc_smooth,row*column,page);
test = permute(R_out_trunc_smooth,[2,3,1]); %reorder the simulated dataset to be dimension x dimension x wavelength
R_out_trunc_smooth_r = reshape(test,size(test,1)*size(test,2),size(test,3));
confidence_r = confidence(:);

t_best_r = zeros(row*column,1);
sig_best_r = t_best_r;
mse_best_r = t_best_r;

test_pad = padarray(test,[2,2],nan,'both');  %create a padded array for performing direct fit. This offsets indices of sigma and thickness by 2 

%This loop takes each measured spectrum from each pixel, replicates the
%spectrum vector into a matrix of similar dimensions as the simulated data
%set, and calculates the MSE of the derivatives between the two
%matrices(i.e., take one spectrum vector, make it m x n to match
%R_out_trunc_smooth_r, take the derivatives of the spectra, and calculate
%the MSE.)
%A row vector of MSEs are calculated. Find the minimum value, and
%find its vector index. This index is used to determine the row and column
%from the simulated data set. The row value corresponds to the matrix index
%of the best thickness value. The column value corresponds to the matrix
%index of the best sigma value. These are then saved in t_best_r,
%sigma_best_r, etc. Outside the loop, these vectors are reshaped into
%matrices with the appropriate dimensions based on the hyperspectral image.
%This method allows for parfor utilization and 50x faster than the original
%script!!
%
tic
parfor i = 1:size(R_meas_trunc_smooth_r,1)
    extract = repmat(R_meas_trunc_smooth_r(i,:),size(R_out_trunc_smooth_r,1),1);
    mse = mean((diff(extract,1,2)-diff(R_out_trunc_smooth_r,1,2)).^2,2);
    %  mse = mean((extract-R_out_trunc_smooth_r).^2,2);
    min_mse = min(mse);%lowest MSE between simulated and measured derivatives
    index_min = find(mse == min(mse));
    if mod(index_min,t_dim)~=0 && floor(index_min/t_dim) ~=0
        index_t = mod(index_min,t_dim);
        index_sig = floor(index_min/t_dim)+1;
    else
        if floor(index_min/t_dim) >0 && mod(index_min,t_dim) ==0
            index_sig = floor(index_min/t_dim);
            index_t = t_dim;
        else
            if floor(index_min/t_dim)==0
                index_sig = 1;
                index_t = mod(index_min,t_dim);
            end
            
        end
    end
    %Second test of direct method +/- 20 nm from t_best and sig_best
    %found from derivative method. This uses a padded array for the
    %simulated data with 2 rows and columns of NaNs to simplify fitting.
%     if index_t+2>= 3 && index_t+2<= t_dim-2 && index_sig+2 >=3 && index_sig+2<=sig_dim-2  %This conditional not necessary
        sub_R_out_trunc_smooth = test_pad(index_t:index_t+4,index_sig:index_sig+4,:);
        sub_R_out_trunc_smooth_r = reshape(sub_R_out_trunc_smooth,size(sub_R_out_trunc_smooth,1)*size(sub_R_out_trunc_smooth,2),size(sub_R_out_trunc_smooth,3));
        extract = repmat(R_meas_trunc_smooth_r(i,:),size(sub_R_out_trunc_smooth_r,1),1);
        mse = mean((extract-sub_R_out_trunc_smooth_r).^2,2);
        index_min = find(mse == min(mse));
        
        if mod(index_min,5)~=0 && floor(index_min/5) ~=0
            index_t_sub = mod(index_min,5);
            index_sig_sub = floor(index_min/5)+1;
        else
            if floor(index_min/5) >0 && mod(index_min,5) ==0
                index_sig_sub = floor(index_min/5);
                index_t_sub = 5;
            else
                if floor(index_min/5)==0
                    index_sig_sub = 1;
                    index_t_sub = mod(index_min,5);
                end
                
            end
        end
       index_t = index_t +(index_t_sub-3);
       index_sig = index_sig+(index_sig_sub-3);
   % end
    t_best_r(i) = t_nacre(index_t);
    sig_best_r(i) = sig_nacre(index_sig);
    mse_best_r(i) = min_mse;
    confidence_r(i) = confidence_r(i)./min_mse;
    best_index{i} = [index_t,index_sig];
end
    toc
t_best = reshape(t_best_r,row,column);
sig_best = reshape(sig_best_r,row,column);
mse_best = reshape(mse_best_r,row,column);
confidence = reshape(confidence_r,size(confidence,1),size(confidence,2));
% confidence =sqrt(confidence); %can turn this back on if needed to scale
% confidence2 = confidence2./mse_best; OPTION: alternate confidence metrics
% confidence3 = 1./mse_best;
best_index = reshape(best_index,row,column);

%Map MTT best fit values
figure                
imagesc(t_best)
axis equal
% colormap(c)
xlabel('x')
ylabel('y')
zlabel('mean T (nm)')
title('Nacre TT Fit')

%Map sigma best fit values
figure
imagesc(sig_best)
axis equal
xlabel('x')
ylabel('y')
zlabel('\sigma (nm)')
title('Nacre \sigma Fit')

%Map best fit MSE
figure
imagesc(mse_best)
axis equal
colormap jet
xlabel('x')
ylabel('y')
zlabel('MSE')
title('MSE Measured v. Fit')

%Map confidence figure of merit
figure
imagesc(confidence)
shading flat
axis equal
colormap jet
xlabel('x')
ylabel('y')
zlabel('Confidence Ratio (Intensity/MSE)')
title('Confidence Ratio')

%% plot best fit spectra and original spectra
% OPTIONAL: Plot the reflectances of the measured and best-fit simulated
% maps to check fitting.

% Click two points to define an area to grab spectra
% for each pixel
% [x,y] = getpts;
% x = round(x);
% y = round(y);
% xy = zeros((abs(x(2)-x(1)))*(abs(y(1)-y(2))),2);
% count = 0;
% for i = 1:abs(x(2)-x(1))
%     for j = 1:abs(y(2)-y(1))
%         count = count+1;
%         xy(count,:) = [x(1)+i-1,y(1)+j-1]; %change to y(1) for imagesc() point grab or y(2) for surface() point grab
%     end
% end
% 
% %Plot individual spectra for each pixel in selected area
% figure
% for i = 1:size(xy,1)
%     temp = R_meas_trunc(xy(i,2),xy(i,1),:);
%                 temp = temp(:); %restructure meas. R data into column vector
%                 
%     temp2 = R_out_trunc(:,best_index{xy(i,2),xy(i,1)}(1),best_index{xy(i,2),xy(i,1)}(2));
%     peak = max(test);
%     test_norm = test./peak; %normalize all reflectance maximums to 1
% %     plot(lambda_trunc(1:end-1), diff(smoothdata(temp2(:),'movmean',10)))  %Plots derivative of simulated spectra that are fit
% %     hold on
% %     plot(lambda_trunc(1:end-1),(diff(smoothdata(temp,'movmean',10))),'k') %Plots derivative of measured spectra that are fit
%     plot(lambda_trunc, temp2(:))
%     hold on
%     plot(lambda_trunc,temp,'k')
% end

%% Data scrubbing
k = 0:0.01:0.99; %scrubbing percentage factor
[xs,ys] = getpts; % select points to define an area that oulines the entire shell. This will be used average points only in shell region, not periphery.

% Use when outlining an area to scrub
t_best_scrub = nan(size(t_best));
sig_best_scrub = nan(size(sig_best));
confidence_scrub = nan(size(confidence)); %For storing only confidence values selected in the region of interest
[x_con,y_con] = find(confidence>=0); %Find x,y indices for all values of original confidence map

%Loop to find only confidence map values within selected region of interest
for i=1:length(x_con)
    if inpolygon(x_con(i),y_con(i),xs,ys)==1
        confidence_scrub(x_con(i),y_con(i)) = confidence(x_con(i),y_con(i));
    end
end

    
% % Use when not defining a specfic area to scrub
% t_best_scrub = t_best;
% sig_best_scrub = sig_best;

% confidence = confidence./sqrt(confidence); %JUST FOR TESTING

for j= 1:length(k)
%     Use when outlining an area to scrub
        conlim =k(j)*max(max(confidence_scrub));  %set the lower bound on the confidence limit
        [x1,y1] = find(confidence_scrub >= conlim);
        t_best_scrub = nan(size(t_best));
        sig_best_scrub = nan(size(sig_best));
        for i = 1:length(x1)
            if inpolygon(x1(i),y1(i),ys,xs) == 1
            t_best_scrub(x1(i),y1(i)) = t_best(x1(i),y1(i));
            sig_best_scrub(x1(i),y1(i)) = sig_best(x1(i),y1(i));
            end
        end
   
%     conlim =k(j)*max(max(confidence));  %set the lower bound on the confidence limit
%     conlim = k(j)*728030;
%     conlim = k(j)*sqrt(133690);
%     [x1,y1] = find(confidence <= conlim);
%     for i = 1:length(x1)
%         t_best_scrub(x1(i),y1(i)) = nan;
%         sig_best_scrub(x1(i),y1(i)) = nan;
%     end
if j==1
    size_unscrub = size(t_best_scrub(~isnan(t_best_scrub))); %Defines the size of unscubbed image to be equal to the reduced area outlined
    end
if j ==2 %|| j == 2 || j == 3 || j == 11  %This if statement just to plot a certain index of scrubbing
    figure
    imagesc(t_best_scrub,'AlphaData',~isnan(t_best_scrub))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('mean T (nm)')
    title('mean T fit')
    caxis([150,600])
    axis off
    colorbar
    title(['percent scrub ',num2str(k(j).*100),'%'])

    figure
    imagesc(sig_best_scrub,'AlphaData',~isnan(t_best_scrub))
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('mean T (nm)')
    title('Scrubbed \sigma Fit')
    colormap jet
    caxis([5,70])
    axis off
    colorbar

 end
    mean_t_scrub(j,1) = nanmean(nanmean(t_best_scrub));
    mean_sig_scrub(j,1) = nanmean(nanmean(sig_best_scrub));
    size_scrub = size(t_best_scrub(~isnan(t_best_scrub)));
    %size_unscrub = size(t_best(:));  %use when not looking at selected
    %area
    image_ratio(j,1) = size_scrub(1)/size_unscrub(1); 
    image_pixels(j,1) = size_scrub(1); %Raw number of pixels used in the sweep average
    
end

figure
plot(k,mean_t_scrub)
figure
plot(k,mean_sig_scrub)
figure
plot(k,image_ratio)
figure
plot(k,image_pixels)
%% Part C: Data analysis

%After performing the scrubbing sweep, place the mean tablet thickness,
%standard deviations, and image ratio for each sample into a matrix. This
%code automatically appends the data to the appropriate matrix.

t_sweep = [t_sweep,mean_t_scrub];
sig_sweep = [sig_sweep,mean_sig_scrub];
img_sweep = [img_sweep,image_ratio];
img_sweep_pixels = [img_sweep_pixels,image_pixels];


