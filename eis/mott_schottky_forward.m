%____________________________________
%
% Author        - Vitamin-C
%
% Status        - Functional
%
% Description   - Models and extracts mott-schottky-like semiconductive
%                 behaviour of examined material based on EIS (electrochemical
%                 impedance spectroscopy results.
% 
% Use Comments  - Set up save dir for first time use
%                 Use find and replace to change variable fields.
%                 This script could be improved with descriptive variable
%                 names.
%
% Initialize parameters and variables:

clear all
clc
v_min=-0.5; %lower limit on potential
v_max=0.3; %upper limit on potential
n_s=((v_max-v_min)/0.025)+1; % number of scans, dependant on potential window
n_s=round(n_s); 
data=zeros(n_s,3); %this is an array of zeros in which to slot our data, the first term is the rows, second is columns, in this case EIS was taken in 33 files, I have 8 variables and wish to record Capacitance as well
Potential=linspace(v_min, v_max, n_s);
m=0.00089196; %scaled sample mass here (usually SSA is used, but mass is simpler and I prefer F/g to F/cm2)
% File - AC        % File prelude
% Folder - EIS\Untreated % Save folder
   
%______________________________________________________________________________
 
for k=1:n_s % loop to create temporary names for each file
    myfilename=sprintf('%s_%d','C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\AC\EIS_1\AC_',k); %set filenames, directory goes here
    A=importdata(myfilename, ' ');%Import as a matrix delimited by the second term (a space for this data)
    A=A.data; %Get rid of those pesky heading lines
    time=A(:,1); %separate and adjust variables
    Zr=A(:,2); %modelling is problematic for small values, so we scale here to milliOhm/g
    Zi=A(:,3);
    Zim=(-1*Zi);
    Hz=A(:,4);
    rads=2*pi*Hz;
    Zcom=sqrt(((Zr).^2)+((Zi).^2)); %This is the magnitude of the complex electric impedance (in Ohm.g)
    Cf=abs((1./(rads.*Zim))); %This is the frequency dependant capacitance of the material (F/g)
    Cf2=Cf.^(-2);
    pote=ones(length(Hz),1)*Potential(1,k);
    tar=strcat('data',num2str(k));
    variable.(tar)=[Hz,pote,Cf2];
end
 
close all
btem=variable.data1;
plot3(btem(:,1),btem(:,2),btem(:,3),'-k')
hold on
 
for h=2:n_s
    nam=strcat('data', num2str(h));
    b=variable.(nam);
    plot3(b(:,1),b(:,2),b(:,3),'-k')
    box on
    xlim([0 21000]);
ylim([-0.55 0.35]);
grid on
view(300,10);
xlabel('Frequency (Hz)')
ylabel('Potential (V vs SCE)')
zlabel('C^{-2} (F^{-2}g^{2})')
set(gca,'fontname','arial')
set(gca,'fontsize',10)
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Nonnormmott\AC(lines).fig');
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Nonnormmott\AC(lines).bmp');
   tem=struct2cell(variable);
   tempy=[tem{1,1};flipud(tem{2,1});tem{3,1};flipud(tem{4,1});tem{5,1};flipud(tem{6,1});tem{7,1};flipud(tem{8,1});...
       tem{9,1};flipud(tem{10,1});tem{11,1};flipud(tem{12,1});tem{13,1};flipud(tem{14,1});tem{15,1};flipud(tem{16,1});tem{17,1};flipud(tem{18,1});tem{19,1};flipud(tem{20,1});...
       tem{21,1};flipud(tem{22,1});tem{23,1};flipud(tem{24,1});tem{25,1};flipud(tem{26,1});tem{27,1};flipud(tem{28,1});tem{29,1};flipud(tem{30,1});tem{31,1};flipud(tem{32,1});tem{33,1}];
end

te=cell2mat(tem);
x=te(:,1);
y=te(:,2);
z=te(:,3);
f = scatteredInterpolant(x,y,z);
xlin = linspace(min(x),max(x),75);
ylin = linspace(min(y),max(y),75);
[X,Y] = meshgrid(xlin,ylin);
Z = f(X,Y);

figure
hold on
surf(X,Y,Z) %interpolated
box on
grid on
xlim([0 21000]);
ylim([-0.55 0.35]);
zlim([0 1.1*max(te(:,3))]);
view(300,10);
xlabel('Frequency (Hz)')
ylabel('Potential (V vs SCE)')
zlabel('C^{-2} (F^{-2}g^{2})')
set(gca,'fontname','arial')
set(gca,'fontsize',10)
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Nonnormmott\AC(surf).fig');
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Nonnormmott\AC(surf).bmp');

sa1=zeros(1,15);
sa2=zeros(15,33);
sa3=zeros(15,33);
namt=strcat('data', num2str(1));
bt=variable.(namt);
sa1=bt(1:15,1);
clear n

for n=1:33
    nam1=strcat('data', num2str(n));
    b1=variable.(nam1);
    
    for k=1:15
    sa3(k,n)=b1(k,3);
    end
    
end

for jl=1:15
    figure
    hold on 
    box on 
    plot(Potential,sa3(jl,:), '-ok')
    xlabel('Potential (V vs SCE)')
    ylabel('C^{-2} (F^{-2}g^{2})')
    set(gca,'fontname','arial')
    set(gca,'fontsize',12)
    Hza=num2str(sa1(jl,1));
    hzn=strcat(Hza,'Hz');
    patha=strcat('C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Nonnormmott\AC(',hzn,').fig');
    pathb=strcat('C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Nonnormmott\AC(',hzn,').bmp');
    saveas(gcf, patha);
    saveas(gcf, pathb);
end

pause(3)
close all