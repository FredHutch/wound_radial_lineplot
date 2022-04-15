function C=wound_radial_lineplot(tiffile,frameid, normtf)
% This function computes 50 radially distributed 1D intensity profiles
% across the wound's center of a syncitial drosophila embryo as described
%in Hui et all 2022 (S. Parkhurst lab). It is assumed that the wound's and
%image's centers coordinates are identical (or very close)
%
% INPUT:            tiffile- path to the image file in .tif format (usually
%                            a time series stack)
%                   frameid- index of frame to analyze
%                   normtf- boolean option for raw intensity (0) or
%                           z-normalized intensity (1)
%
%OUTPUT:            C- a 301-by-50 array, one intensity profile per column,
%                      one pixel per row, pixel 151 being the center of the
%                      wound/image
%
% Results are also plotted: x-axis, distance from center in pixels; y-axis,
% fluorescent intensity; thin gray lines, individual profiles; thick black
% line, mean profile.
%
% Example:
%
% C=wound_radial_lineplot('WoundExample.tif',8, 1);
%
% Julien Dubrulle, 2022

%% Gather image info, read file and smooth data

close all;
INFO=imfinfo(tiffile);
if frameid>numel(INFO)
    clc;
    disp(['FRAME OUT OF RANGE (' int2str(numel(INFO)) ' TOTAL)'])
end
x=INFO.Width;
y=INFO.Height;
I=imread(tiffile,frameid);
ISg=imgaussfilt(I,2);

%% Initialize parameters
diami=linspace(0,300,301); % size of the diameter of circle around center of wound (center of image)
cent=[x/2 y/2]; % center of image/Wound
da=linspace(0,pi,51); %angles defining all orientations analyzed (def=50)
dan=-(da); %get opposite/antipodal points
[xc, yc]=pol2cart([da fliplr(dan)],repmat(150,1,102)); %get cartesian subscripts of points with radius 150
xx=xc+cent(1); %use center of image as ref
yy=yc+cent(2);

C=zeros(301,51); %Intialize traces)


%% get trace for each diameter/line
for i=1:51
    Ct=improfile(ISg,[xx(i) xx(i+51)],[yy(i) yy(i+51)]); %get trace
    
    %Make sure traces have the same number of elements- interpolate if
    %necessary
    if numel(Ct)<301
        diam=linspace(0,300,numel(Ct));
        Cf=interp1(diam,Ct,diami);
        C(:,i)=Cf;
    else
        C(:,i)=Ct;
    end
end

C(:,end)=[]; %redondant profile removed

%% Z-Normalization
if normtf
    C=(C-mean(C))./std(C);
end

%% Plot result
diamix=diami-151;

figure,plot(diamix,C,'color',[0.5 0.5 0.5])
hold on;
Cm=mean(C,2);
plot(gca,diamix,Cm,'-k','Linewidth',2);
if normtf
    ylabel('z-norm. signal intensity')
else
    ylabel('Signal intensity');
end
xlabel('Distance from wound center (in pixels)');
xlim([-150 150]);

