%____________________________________

% Author        - Vitamin-C
% Status        - Functional
% Description   - Process cathodic and anodic cycling
%                 over n cycles for a set potential range. Plots the
%                 voltammogram, the absolute variation in integrative
%                 capacitance with cycle number, and the fractional
%                 variation in integrative capacitance with cycle number.
% Use Comments  - Set up save dir for first time use
%                 Files being read are without extention due to the manner
%                 exported by iviumstat device, likely requires alteration
%                 for any other files.
%                 Use find and replace to change variable fields.
%                 This script could be improved with descriptive variable
%                 names.
% Initial Prep
clear all
clc
close all
%____________________________________

% Set the variable fields for your sample

% Save Directory - CV retested\Copper
% Filename - Cu_1_org_2
% Mass - 0.00056705 % Electroactive mass of the sample.
n = 1000; % Total number of cycles.
v_max = 1; % Upper potential limit
v_min = -1.5; % Lower potential limit
m = 950; % How many end cycles calculate the capacitance.

%____________________________________

h = figure;
for k = 20 : n %n here

% Create sequential filenames from base filename
myfilename = sprintf('%s_%d', 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\Cu_1_org_2\CV_1\Cu_1_org_2_', k);

% Import files
A = importdata(myfilename, ' ', 1);  
B = A.data;

%Assign column vectors
t = B(:, 1);
E = B(:, 2);
i = (B(:, 3)) / 0.00056705; %mass here

%plot potential vs current

%find max current limits
i_max = max(abs(i));
currentmax(k,:) = i_max;

plot(E, i, 'k', 'LineWidth', 0.25)
hold on
end

%find max current limits
box on
i_max = max(abs(currentmax));
i_min = -1 * i_max;
axis([v_min v_max i_min i_max]) %vmin and vmax
set(gca, 'xTick', v_min:0.1:v_max, 'Fontname', 'Arial', 'Fontsize', 12) %vmin and vmax
xlabel('Potential (V vs SCE)','Fontsize',14,'Fontname','Arial')
ylabel('Current (A/g)','Fontsize',14,'Fontname','Arial')

%save plot to specified directory
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CV_1.fig')
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CV_1.bmp')

%create empty array for data to fill
data=zeros(n,3);

for k=1:n 

%Create sequential filenames from base filename
myfilename=sprintf('%s_%d','D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\Cu_1_org_2\CV_1\Cu_1_org_2_',k); %set filename

%import files
B=importdata(myfilename,' ',1);
A=B.data;

%Assign column vectors
t=A(:,1);
E=A(:,2);
i=A(:,3);
 
%Extract data points where i>0 into Anodic matrix  
%Extract data points where i<0 into Cathodic matrix
Anodic=A(A(:,3)>0,:);
Cathodic=A(A(:,3)<0,:);
t_a=Anodic(:,1);
E_a=Anodic(:,2);
i_a=Anodic(:,3);

%integrate anodic current (i_a) wrt t
Z=cumtrapz(t_a,i_a);
Capacity_Anodic=((Z(end)/0.00056705))/(v_max-v_min); %mass, v_max, v_min here
 
%integrate cathodic current (i_c) wrt t
t_c=Cathodic(:,1);
E_c=Cathodic(:,2);
i_c=Cathodic(:,3);
Y=cumtrapz(t_c,i_c);
Capacity_Cathodic=(abs(Y(end))/0.00056705)/(v_max-v_min); %mass, v_max, v_min here
 
%add results to data matrix
data(k,:)=[k;Capacity_Anodic;Capacity_Cathodic];
end

%show results array
Results=data;
Cycle=Results(:,1);
C_A=Results(:,2);
C_C=Results(:,3);
max_An=max(C_A);
max_Ca=max(C_C);
max_cap=max(max_An,max_Ca);
 
last_m=n-m; 
 
A_last=C_A(last_m:n); 
C_last=C_C(last_m:n); 
 
Av_A=mean(A_last);
Av_C=mean(C_last);
AV=(Av_A+Av_C)/2;
  
%Plot results 
r=figure;
hold on
box on
plot(Cycle,C_A,'k','LineWidth',2)
axis([0 n 0 250])
set(gca,'xTick',0:(n/10):n,'Fontsize',12,'Fontname','Arial') %n here
set(gca,'YLim',[0 250])

hold on % Pointless line, confirm later before removing

plot(Cycle,C_C,'Linestyle','--','Color','k','Linewidth',2)
hold off
str1=sprintf('Anodic [%0.1f] (F/g)',Av_A);
str2=sprintf('Cathodic [%0.1f] (F/g)',Av_C);
rleg=legend(str1,str2);
xlabel('Cycle number','Fontsize',14,'Fontname','Arial')
ylabel('Capacitance (F/g)','Fontsize',14,'Fontname','Arial')
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CS_1.fig');
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CS_1.bmp')

%Plot fractional results 
h=figure;
hold on
box on
plot(Cycle,(C_A./C_A(last_m))*100,'k','LineWidth',2)
axis([0 n 0 250]) %n here and below
set(gca,'xTick',0:(n/10):n,'Fontsize',12,'Fontname','Arial') %n here
set(gca,'YLim',[0 250])
hold on
plot(Cycle,(C_C./C_C(last_m))*100,'Linestyle','--','Color','k','Linewidth',2)
hold off
str1=sprintf('Anodic [%0.1f] (F/g)',Av_A);
str2=sprintf('Cathodic [%0.1f] (F/g)',Av_C);
rleg=legend(str1,str2);
xlabel('Cycle number','Fontsize',14,'Fontname','Arial')
ylabel('Fractional Capacitance (%)','Fontsize',14,'Fontname','Arial')
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CS%_1.fig');
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CS%_1.bmp')

s=figure;
for k=20:200 %I have 200 cycles in my second set
      %Create sequential filenames from base filename
myfilename=sprintf('%s_%d','D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\Cu_1_org_2\CV_2\Cu_1_org_2_',k); %filename here
%import files
A=importdata(myfilename,' ',1);
B=A.data;
%Assign column vectors
t=B(:,1);
E=B(:,2);
i=(B(:,3))/0.00056705; %mass here
%plot E vs i
%find max current limits
i_max=max(abs(i));
currentmax(k,:)=i_max;

plot(E,i,'k','LineWidth',0.25)
hold on
end
%find max current limits
box on
i_max=max(abs(currentmax));
i_min=-1*i_max;
axis([v_min v_max i_min  i_max]) %vmin and vmax
set(gca,'xTick',v_min:0.1:v_max,'Fontname','Arial','Fontsize',12)%vmin and vmax
xlabel('Potential (V vs SCE)','Fontsize',14,'Fontname','Arial')
ylabel('Current (A/g)','Fontsize',14,'Fontname','Arial')
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CV_2.fig')
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CV_2.bmp')
%save plot to specified directory
 
%Manual save to directory
%MATLAB code used to calculate area under CV curve and to plot capacitance against sweep number.

%create empty array for data
data=zeros(200,3); %n here
for k=1:200 %n here
    %Create sequential filenames from base filename
myfilename=sprintf('%s_%d','D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\Cu_1_org_2\CV_2\Cu_1_org_2_',k); %set filename
%import files
B=importdata(myfilename,' ',1);
A=B.data;
%Assign column vectors
t=A(:,1);
E=A(:,2);
i=A(:,3);
 
%Extract data points where i>0 into Anodic matrix  
%Extract data points where i<0 into Cathodic matrix
Anodic=A(A(:,3)>0,:);
Cathodic=A(A(:,3)<0,:);
t_a=Anodic(:,1);
E_a=Anodic(:,2);
i_a=Anodic(:,3);
%integrate anodic current (i_a) wrt t
Z=cumtrapz(t_a,i_a);
Capacity_Anodic=((Z(end)/0.00056705))/(v_max-v_min); %mass, v_max, v_min here
 
%integrate cathodic current (i_c) wrt t
t_c=Cathodic(:,1);
E_c=Cathodic(:,2);
i_c=Cathodic(:,3);
Y=cumtrapz(t_c,i_c);
Capacity_Cathodic=(abs(Y(end))/0.00056705)/(v_max-v_min); %mass, v_max, v_min here
 
%add results to data matrix
data(k,:)=[k;Capacity_Anodic;Capacity_Cathodic];
end
%show results array
Results=data;
Cycle=Results(:,1);
C_A=Results(:,2);
C_C=Results(:,3);
max_An=max(C_A);
max_Ca=max(C_C);
max_cap=max(max_An,max_Ca);
 
last_m=200-150; %n here
 
A_last=C_A(last_m:200); %n here
C_last=C_C(last_m:200); %n here
 
Av_A=mean(A_last);
Av_C=mean(C_last);
AV=(Av_A+Av_C)/2;
  
%Plot results 
q=figure;
hold on
box on
plot(Cycle,C_A,'k','LineWidth',2)
axis([0 200 0 250]) %n here and below
set(gca,'xTick',0:(200/10):200,'Fontsize',12,'Fontname','Arial') %n here
%set(gca,'YLim',[0 (max_cap*1.3)])
set(gca,'YLim',[0 250])
hold on
plot(Cycle,C_C,'Linestyle','--','Color','k','Linewidth',2)
hold off
str1=sprintf('Anodic [%0.1f] (F/g)',Av_A);
str2=sprintf('Cathodic [%0.1f] (F/g)',Av_C);
qleg=legend(str1,str2);
xlabel('Cycle number','Fontsize',14,'Fontname','Arial')
ylabel('Capacitance (F/g)','Fontsize',14,'Fontname','Arial')
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CS_2.fig');
saveas(gcf, 'D:\Documents\OneDrive\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Copper\Cu_1_org_2_CS_2.bmp')
%save plot to specified directory
close all
%Manual save to directory

