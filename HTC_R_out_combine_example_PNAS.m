clear all;
close all;

%Combine HTC generated output files for the R_out files for
%501 randomizations, p-polarized, 30 degree incident angle,with high
%magnification (i.e., +-9 degrees incident angle). Jad Salman.
path = pwd;

%Create a structure with all output file names from fracture HTC output
%files
output = dir(fullfile(path,'out*.*')); 

%Initialize a matrix of the same dimensions of the desired R_out file used
%to fit data in nacre model
total_out = zeros(204,91,14); 

%Loop loads every output file name incrementally into 'temp', then finds
%the components of the output file name that are numeric (e.g. 10 in
%'out10.mat') and increments it by one to match the true column number the
%set of data belongs to. The main point is to ensure that the structure of
%mean thickness (91 columns) and standard deviation (14 pages) is in the
%proper order for running the nacre fitting script.

for i = 1:91
temp = load(output(i).name);
% output(i).name
column = str2double(output(i).name(4:5))+1; %finds proper column number
total_out(:,column,:) = temp.R_out;
end

R_out = total_out;
