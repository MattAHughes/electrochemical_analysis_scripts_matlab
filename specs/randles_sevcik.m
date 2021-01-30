%____________________________________

% Author        - Vitamin-C
% Status        - Functional
% Description   - Process step potential electrochemical spectroscopy
%                 (SPECS) results similar to cycling rate then apply the
%                 randles sevcik equation and calculate surface coverage by
%                 hydrogen ions.
% Use Comments  - Set up save dir for first time use
%                 Use find and replace to change variable fields.
%                 This script could be improved with descriptive variable
%                 names.

% Initial Prep
clc
clear all
close all

% Finding and replace a systematic filename

% filename - 700_4
% mass - 0.00038519 % Electroactive mass
% directory - SPECS\Untreated % Where files are stored
% saves - SPECS Results\Adsorbed % Where output will be saved
%___________________________________________________________________

cd('C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS\Untreated')

filenamea='700_4_1.csv'; %Forward Sweep
delimiterIn=',';
A=importdata(filenamea,delimiterIn);
filenameb='700_4_2.csv'; %Back Sweep
C=importdata(filenameb,delimiterIn);

%______________________________________

E=0.025; %Potential Step
srate=(E)./linspace(0.2,299.8,1499); %scan rate in V s^-1
time=299.8; %How long you wish to look at
mass=0.00038519; %parameters go here
Czero=0.5; %original conc, mol/L
Fara=96485;
nele=1; 

%______________________________________

T1=A(:,1); % to plot the total SPECS experiment separate time and current
I1=A(:,2);
i1=I1; 
T2=C(:,1); 
I2=C(:,2);
i2=I2;
t2=T2+10200; %shift across so both can be plotted together
Potential=linspace(-0.500,0.300, 33)';
Potentialr=linspace(0.275,-0.475,31);
x=linspace(0,20299,20300);
y=-0.525+0.025*(heaviside(x-300)+heaviside(x-600)+heaviside(x-900)+heaviside(x-1200)+heaviside(x-1500)+heaviside(x-1800)+heaviside(x-2100)+heaviside(x-2400)+heaviside(x-2700)+heaviside(x-3000)+heaviside(x-3300)+heaviside(x-3600)+heaviside(x-3900)+heaviside(x-4200)+heaviside(x-4500)+heaviside(x-4800)+heaviside(x-5100)+heaviside(x-5400)+heaviside(x-5700)+heaviside(x-6000)+heaviside(x-6300)+heaviside(x-6600)+heaviside(x-6900)+heaviside(x-7200)+heaviside(x-7500)+heaviside(x-7800)+heaviside(x-8100)+heaviside(x-8400)+heaviside(x-8700)+heaviside(x-9000)+heaviside(x-9300)+heaviside(x-9600)+heaviside(x-9900)-(heaviside(x-10200)+heaviside(x-10500)+heaviside(x-10800)+heaviside(x-11100)+heaviside(x-11400)+heaviside(x-11700)+heaviside(x-12000)+heaviside(x-12300)+heaviside(x-12600)+heaviside(x-12900)+heaviside(x-13200)+heaviside(x-13500)+heaviside(x-13800)+heaviside(x-14100)+heaviside(x-14400)+heaviside(x-14700)+heaviside(x-15000)+heaviside(x-15300)+heaviside(x-15600)+heaviside(x-15900)+heaviside(x-16200)+heaviside(x-16500)+heaviside(x-16800)+heaviside(x-17100)+heaviside(x-17400)+heaviside(x-17700)+heaviside(x-18000)+heaviside(x-18300)+heaviside(x-18600)+heaviside(x-18900)+heaviside(x-19200)+heaviside(x-19500)+heaviside(x-19800)+heaviside(x-20400)));
sca=268600*(1.84095e-5^0.5)*0.0005;

%___________________________________ correct for double-stacked excitations


%Next we examine curves in the forward sweep

%_____________________________________________________________________________________

it=zeros(length(srate),33); % storage matrix
for k=1:33                        %Forward Sweep
    B=A(((1500*k)+1):(1500*k+time*5),:);
    t=(linspace(0.2, time, time*5))';
    i=i1((1500*k+1):(1500*k+time*5),:);
   %take i at each t
    for j1=1:length(i)
        it(j1,k)=(i(j1)*4*8.314*298.15)/(96485.332^2);
    
    end
  
end

itp=[Potential';it];
itr=zeros(length(srate),31);

for z=1:31                        
    D=C(1500*z+1:1500*z+time*5,:);
    t=(linspace(0.2, time, time*5))';
    i=1*i2(1500*z+1:1500*z+time*5,:);
    for j2=1:length(i)
        itr(j2,z)=(i(j2)*4*8.314*298.15)/(96485.332^2);
    
    end
end

figure
hold on
box on

for y=[1 4 7 10 40 70 100 400 700 1000 1400]
    plot([Potential' Potentialr Potential(1,1)],smooth([it(y,:) itr(y,:) it(y,1)],'rlowess'),'k')
    txt=sprintf('%f V s^{-1}',srate(1,y));
    text(-0.2,max(smooth([it(y,:) itr(y,:) it(y,1)],'rlowess')),txt,'fontname','arial','fontsize',6,'backgroundcolor','w');
  
end

xlabel('Potential (V vs SCE)')
ylabel('{4IRT}/{n^2F^2} (mol V s^{-1})')
set(gca,'fontname','arial')
set(gca,'fontsize',12)
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Adsorbed\700_4_CVscanrate.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Adsorbed\700_4_CVscanrate.bmp');
sr=srate;
xx=linspace(sr(1),sr(end),5000);
ft=fittype('linear');
m=zeros(33,4);
b=zeros(33,4);
figure
plot(srate,it,'-k')
xlabel('Scan Rate (V s^{-1})')
ylabel('{4IRT}/{n^2F^2} (mol V s^{-1})')
set(gca,'fontname','arial')
set(gca,'fontsize',12)
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Adsorbed\700_4_ivsrate.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Adsorbed\700_4_ivsrate.bmp');

for n=1:33
    var={sprintf('%f V',Potential(n,1))};
  cf=fit(sr',it(:,n),ft);
  fita=fitlm(xx,cf(xx));
  m(n,:)=table2array(fita.Coefficients(2,:));
  b(n,:)=table2array(fita.Coefficients(1,:));
figure
plot(fita)
text(min(xx),max(it(:,n)),var,'fontname','arial','fontsize',6,'backgroundcolor','w');
xlabel('scan rate (V s^{-1})')
ylabel('{4IRT}/{n^2F^2} (mol V s^{-1})')
set(gca,'fontname','arial')
set(gca,'fontsize',12)
end

header={'Slope','SE','tstat','pvalue'};
header2={'Intercept','SE','tstat','pvalue'};
T = table(m(:,1),m(:,2),m(:,3),m(:,4));
T.Properties.VariableNames = header;
T2 = table(b(:,1),b(:,2),b(:,3),b(:,4));
T2.Properties.VariableNames = header2;
plot(Potential,smooth(m(:,1)./mass,'rlowess'),'-k')
hold on
plot(Potential,m(:,1)./mass,'-.k')
xlabel('Potential (V vs SCE)')
ylabel('Surface Coverage (mol g^{-1})')
set(gca,'fontname','arial')
set(gca,'fontsize',12)
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Adsorbed\700_4_adsorbedmolg-1.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Adsorbed\700_4_adsorbedmolg-1.bmp');

pause(3)
close all
