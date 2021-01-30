%____________________________________

% Author        - Vitamin-C
% Status        - Functional
% Description   - Process cathodic and anodic cycling
%                 over n cycles for a symmetrically increasing potential range. Plots the
%                 voltammogram and the absolute variation in integrative
%                 capacitance over a set number of cycles at each potential step.
% Use Comments  - Set up save dir for first time use
%                 Loops over folders in a directory, each of which contains
%                 several individual (extentionless) cv files (iviumstat
%                 allows us to export with no extention).
%                 Use find and replace to change variable fields.
%                 This script could be improved with descriptive variable
%                 names.
% Initial Prep
clear all
clc
%____________________________________

% Save Directory - CV retested\Doped
% Filename - nitrate_na2so4_W
% Mass - 0.00054184
n=5; % number of cycles
m=0; %how many end cycles calculate capacitance
cv_n=33; %the number of folders with CVs in them

h=figure;
for l=1:cv_n
dir=sprintf('%s_%d','C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\nitrate_na2so4_W\CV',l);
name=strcat(dir,'\nitrate_na2so4_W_');

for k=2:n %n here
%Create sequential filenames from base filename
myfilename=sprintf('%s_%d',name,k); %filename here

%import files
A=importdata(myfilename,' ',1);
B=A.data;

%Assign column vectors
t=B(:,1);
E=B(:,2);
i=(B(:,3))/0.00054184; %mass here

%plot E vs i
%find max current limits
i_max=max(abs(i));
currentmax(k,:)=i_max;
v_max=max(E);
v_min=min(E);
plot(E,i,'k','LineWidth',0.05)
hold on
end

end

box on
i_max=max(abs(currentmax));
i_min=-1*i_max;
axis([v_min v_max i_min  i_max]) %vmin and vmax
set(gca,'xTick',v_min:0.1:v_max,'Fontname','Arial','Fontsize',12)%vmin and vmax
xlabel('Potential (V vs SCE)','Fontsize',14,'Fontname','Arial')
ylabel('Current (A/g)','Fontsize',14,'Fontname','Arial')
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Doped\nitrate_na2so4_W_windowCV_1.fig')
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Doped\nitrate_na2so4_W_windowCV_1.bmp')

%create empty array for data
data=zeros(n,3); %n here

r=figure;

for l=1:cv_n   
dir=sprintf('%s_%d','C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\nitrate_na2so4_W\CV',l);
name=strcat(dir,'\nitrate_na2so4_W_');

for k=1:n %n here

%Create sequential filenames from base filename
myfilename=sprintf('%s_%d',name,k); %set filename

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
Capacity_Anodic=((Z(end)/0.00054184))/(v_max-v_min); %mass, v_max, v_min here
 
%integrate cathodic current (i_c) wrt t
t_c=Cathodic(:,1);
E_c=Cathodic(:,2);
i_c=Cathodic(:,3);
Y=cumtrapz(t_c,i_c);
Capacity_Cathodic=(abs(Y(end))/0.00054184)/(v_max-v_min); %mass, v_max, v_min here
 
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
 
last_m=n-m; %n and m here
 
A_last=C_A(last_m:n); %n here
C_last=C_C(last_m:n); %n here
 
Av_A=mean(A_last);
Av_C=mean(C_last);
AV=(Av_A+Av_C)/2;
  
%Plot results 
hold on
box on
plot(Cycle,C_A,'k','LineWidth',2)
axis([0 n 0 250]) %n here and below
set(gca,'xTick',0:(n/5):n+1,'Fontsize',12,'Fontname','Arial') %n here
xlim([0 6])

%set(gca,'YLim',[0 (max_cap*1.3)])
set(gca,'YLim',[0 250])
xlabel('Cycle number','Fontsize',14,'Fontname','Arial')
ylabel('Capacitance (F/g)','Fontsize',14,'Fontname','Arial')
end

saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Doped\nitrate_na2so4_W_windowCS_1.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\CV retested\Doped\nitrate_na2so4_W_windowCS_1.bmp')

pause(5)
close all