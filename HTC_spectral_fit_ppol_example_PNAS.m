function HTC_spectral_fit_ppol_501_rand_himag_14deg(varargin)
% This is a sample code for performing transfer-matrix method + Monte Carlo
% method calculations for simulated nacre stacks. The code is presented in
% a function to allow for high throughput computing (HTC). Jad Salman.
clear all

path = pwd;

%% Define parameters
%Wavelength define:
load('test_reflectances_iq.mat') %sample reflectance data from Specim IQ
lambda=test_reflectances_iq(:,1)*1e-9; %Wavelength [m]
test_R = test_reflectances_iq(:,2); %For comparing to a measured spectrum

%Angle of incidence + angles to average spectra over. Defined by the
%acceptance angle of the optics in your system.
theta= linspace(5,23,5); % 30 degrees +-9degrees incident angle (based on 8 cm distance from 1 inch diameter lense aperture.)

% Load random number generator saved state to seed same values
load('random_seed.mat')
rng(random_seed);

% Load input of thickness for specific job. This is specifc to running on
% HTC. The 'in.*' file is a text file with one row vector of an m x n
% matrix that defines the thicknesses and standard deviations for which to
% calculate spectra. In this case, we have 91 'in' files each with 14
% standard deviations.
input = dir(fullfile(path,'in*.*'));
load(input.name);

%Refractive index of material layers
air =1;
n_aragonite = 1.63; %
n_organic = 1.43;
n_substrate = 1.63; %substrate is semi infinite. Assumed to be aragonite here

%% Define layer thicknesses  in nm
layers = 201; %Number of layers
number_rand = 501; %number of randomly generated layer stacks.
avg_thick_aragonite = in(1)*1e-9; % Set an average target tickness of aragonite. in(1) is the thickness for one block of standard deviations to loop through
avg_thick_organic = 25*1e-9; % Set an average target thickness of organic
tol_aragonite = in(2:end)*1e-9; % Set tolerances [m]
tol_organic = 10e-9; % fixed value for the organic layers [m]

%Initialize reflectance spectra output for each loop

R_out = zeros(length(lambda),length(avg_thick_aragonite),length(tol_aragonite));
% mse_out = zeros(length(avg_thick_aragonite),length(tol_aragonite));
for t = 1:length(avg_thick_aragonite)
    for h = 1:length(tol_aragonite)
% layer_1_t_min = avg_thick_organic-tol_organic; %Layer 1 thickness [m]
% layer_1_t_max = avg_thick_organic+tol_organic;
% layer_2_t_min = avg_thick_aragonite-tol_aragonite; %Layer 2 thickness [m]
% layer_2_t_max = avg_thick_aragonite+tol_aragonite;
thickcell = zeros(layers,number_rand);
%Creates an array of thicknesses while also randomly generating layer stack
%variability within a range of thicknesses
for i=1:number_rand
    for j = 1:layers
        if mod(j,2)==1
            %Take absolute value of randomized thickness to ensure no
            %negative values are generated. Prevents gain in TM method
            %calculation
            thickcell(j,i) = abs(normrnd(avg_thick_aragonite(t),tol_aragonite(h)));%(layer_2_t_max - layer_2_t_min).*rand+layer_2_t_min;
        else
            thickcell(j,i) = abs(normrnd(avg_thick_organic,tol_organic)); %(layer_1_t_max - layer_1_t_min).*rand+layer_1_t_min;
        end
    end
end

% Define layer index matrix
n_layers = zeros(size(thickcell));
for i = 1:layers+2
    if i == 1 %First layer must be air
        n_layers(i,:) = air;
    else if mod(i,2) == 0 && i ~=layers+2
            n_layers(i,:) = n_aragonite;
        else if mod(i,2) == 1 && i ~=layers+2
                n_layers(i,:) = n_organic;
            else
                n_layers(i,:) = n_aragonite;%Last layer must be aragonite
            end
        end
    end
end
%% RT transfer matrix calculation

%Initialize a matrix to contain the entire output of R and T for all
%angles and all randomly generated layers across all wavelengths.
Rp = zeros(length(lambda),number_rand, length(theta));
% Rs = zeros(length(lambda),number_rand, length(theta));
Out_Rp_avg = zeros(length(lambda),length(theta));
% Out_Rs_avg = zeros(length(lambda),length(theta));

%Calculate R and T. Input result into appropriate row,column,page of Out_R
%and Out_T
% tic
for j = 1:number_rand
    
    [Rp(:,j,:)] = trnsfr_mm_dspless(layers+2,thickcell(:,j),n_layers(:,j),lambda,theta);
    %         R=(Out0(:,1)*100+Out0(:,3)*100)/2;
    %         T=(Out0(:,2)*100+Out0(:,4)*100)/2;
    %         Out_R(:,j,i)=R;
    %         Out_T(:,j,i)=T;
    
end
% time(t) = toc

%% Average reflection output at each defined angle over all lambda
% figure(1)
% figure(2)
for i = 1:length(theta)
   Out_Rp_avg(:,i) = sum(Rp(:,:,i),2)./number_rand;
%    Out_Rs_avg(:,i) = sum(Rs(:,:,i),2)./number_rand;
%    figure(1)
%    plot(lambda,Out_Rp_avg(:,i));
%    hold on
%    figure(2)
%    plot(lambda,Out_Rs_avg(:,i));
%    hold on
end
% figure(1)
% title('p-pol')
% figure(2)
% title('s-pol')
%% Average over angles of incidence

%Average over +/- 10 degrees of incidence angle to replicate field of view
% Rp_0_out = sum(Out_Rp_avg(:,1:3),2)./3;
% Rp_20_out = sum(Out_Rp_avg(:,1:4),2)./5;
% Rp_30_out = sum(Out_Rp_avg(:,3:7),2)./5;
% Rp_40_out = sum(Out_Rp_avg(:,5:9),2)./5;
% Rp_50_out = sum(Out_Rp_avg(:,7:11),2)./5;
% Rp_60_out = sum(Out_Rp_avg(:,9:13),2)./5;

Rp_62_out = sum(Out_Rp_avg(:,1:5),2)./5;
% Rs_62_out = sum(Out_Rs_avg(:,1:5),2)./5;

% figure
% hold on
% plot(lambda,Rp_20_out);
% plot(lambda,Rp_30_out);
% plot(lambda,Rp_40_out);
% plot(lambda,Rp_50_out);
% plot(lambda,Rp_60_out);

% plot(lambda,Rp_62_out);
% title('p-pol')
% figure
% plot(lambda,Rs_62_out);
% title('s-pol')
Rp_62_smooth = smooth(Rp_62_out,1); % Smoothing optional, set to 1 here.
Norm_Rp_62_smooth = Rp_62_smooth./max(Rp_62_smooth);
% Norm_test_62 = test_R./max(test_R);

%Mean squared error
% mse_out(t,h) = immse(Norm_test_62,Norm_Rp_62_smooth)

%Save normalized reflection spectra for each loop
R_out(:,t,h) = Norm_Rp_62_smooth;
%% plot of measured and calculated data
% figure(10)
% plot(lambda,Norm_test_62,'k','linewidth',3)
% hold on
% plot(lambda,Norm_Rp_62_smooth)
% hold on
    end
end
save(strcat('out_',input.name(3:end)),'R_out')


