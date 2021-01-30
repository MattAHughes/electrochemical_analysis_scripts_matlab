%____________________________________

% Author        - Vitamin-C

% Status        - Functional

% Description   - Process step potential electrochemical spectroscopy data
%                 produce information on: - Current and potential with step
%                 - Individual step fittings and their residuals - A
%                 breakdown of the current during a step attributed to
%                 double layer, diffusion, and residual effects - An
%                 overview of the total resistive and capacitive effects
%                 with time - A Ragone diagram for the system.
% 
% Use Comments  - Set up save dir for first time use
%                 Use find and replace to change variable fields.
%                 This script could be improved with descriptive variable
%                 names.

% Initial Prep
% filename - Cu_1_org_2
% mass - 0.00056705 % Electroactive mass
% directory - Supercapacitor\Capacitance stuff % Load directory
% saves - SPECS Results\Copper % Save directory
clc
clear all
v_min=-1.5;
v_max=1;

%___________________________________________________________________

filenamea=sprintf('%s_%d','C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\Cu_1_org_2\SPECS_1\Cu_1_org_2',1); %filename here
delimiterIn=' ';
A=importdata(filenamea,delimiterIn);
filenameb=sprintf('%s_%d','C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\Cu_1_org_2\SPECS_2\Cu_1_org_2',1); %filename here; %Back Sweep
C=importdata(filenameb,delimiterIn);
%______________________________________

E=0.025; %Potential Step
time=50; %How long you wish to look at
mass=0.00056705; %parameters go here
Czero=0.5; %original conc, mol/L
Fara=96485;
nele=1; 

%______________________________________

T1=A(:,1); % to plot the total SPECS experiment separate time and current
I1=A(:,2);
i1=I1/mass; 
p1=A(:,3);
T2=C(:,1); 
fna=((max(T1)*0.1)/30)-1; %number of peaks in forward sweep
rn=fna-1; %number of peaks in reverse sweep
I2=C(:,2);
i2=I2/mass;
t2=T2+T1(end); %shift across so both can be plotted together
p2=C(:,3);
y=[p1' p2'];
Potential=linspace(v_min,v_max, 60)';
x=linspace(0,36300,24291);
Potentialtrue=-0.525+0.025*(heaviside(x-300)+heaviside(x-600)+heaviside(x-900)+heaviside(x-1200)+heaviside(x-1500)+heaviside(x-1800)+heaviside(x-2100)+heaviside(x-2400)+heaviside(x-2700)+heaviside(x-3000)+heaviside(x-3300)+heaviside(x-3600)+heaviside(x-3900)+heaviside(x-4200)+heaviside(x-4500)+heaviside(x-4800)+heaviside(x-5100)+heaviside(x-5400)+heaviside(x-5700)+heaviside(x-6000)+heaviside(x-6300)+heaviside(x-6600)+heaviside(x-6900)+heaviside(x-7200)+heaviside(x-7500)+heaviside(x-7800)+heaviside(x-8100)+heaviside(x-8400)+heaviside(x-8700)+heaviside(x-9000)+heaviside(x-9300)+heaviside(x-9600)+heaviside(x-9900)-(heaviside(x-10200)+heaviside(x-10500)+heaviside(x-10800)+heaviside(x-11100)+heaviside(x-11400)+heaviside(x-11700)+heaviside(x-12000)+heaviside(x-12300)+heaviside(x-12600)+heaviside(x-12900)+heaviside(x-13200)+heaviside(x-13500)+heaviside(x-13800)+heaviside(x-14100)+heaviside(x-14400)+heaviside(x-14700)+heaviside(x-15000)+heaviside(x-15300)+heaviside(x-15600)+heaviside(x-15900)+heaviside(x-16200)+heaviside(x-16500)+heaviside(x-16800)+heaviside(x-17100)+heaviside(x-17400)+heaviside(x-17700)+heaviside(x-18000)+heaviside(x-18300)+heaviside(x-18600)+heaviside(x-18900)+heaviside(x-19200)+heaviside(x-19500)+heaviside(x-19800)+heaviside(x-20400)));
%___________________________________ correct for double-stacked excitations

figure
hold on
yyaxis left %start with the left axis plot (requires matlab 2016 edition or later for yyaxis command to work)
plot(T1, i1, 'k', t2, i2, 'k-');
xlim([0 t2(end)]);
ylim([-10 10]);
no=xlabel('Time (s)');
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
set(no,'Fontname', 'Arial', 'Fontsize', 14);
ye=ylabel('Current (A/g)');
set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
set(gca, 'YColor', 'k') %set colour of axis to black
box on

yyaxis right %onwards to the right axis plot (Potential)
plot([T1' t2'],y,'-r');
ze=ylabel('Potential (V vs SCE)');
set(gca, 'YColor', 'k') %set colour of axis to black
set(ze, 'Fontname', 'Arial', 'Fontsize', 14)
hold off

saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_compilation.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_compilation.bmp');

%Next we examine curves in the forward sweep
%_____________________________________________________________________________________

fwd=zeros(fna,7);                 %For variable storage

for k=1:fna                        %Forward Sweep
    
    B=A(((3000*k)+1):(3000*k+time*5),:);
    t=(linspace(0.1, time, time*5))';
    i=i1((3000*k+1):(3000*k+time*5),:);
    lb=[0 0 0 0 0 -5]; % lower bound on constants (whatever is reasonable)
    ub=[10 400 10 300 20 5]; % upper bound on constants
    c0=[3.3005e-02   9.4509e+01   6.5139e-03   8.2222e+01   7.4400e-01  -2.0037e-01]; %Initial estimates (as close as possible to true values, use EIS results etc. if necessary, can get some very unlikely results otherwise)
    
    f1=@(c,t) (E/c(1))*exp(-t./(c(1)*c(2)))+(E/c(3))*exp(-t./(c(3)*c(4)))+(c(5)./(t.^0.5))+c(6);
    options = optimoptions('lsqcurvefit','FunctionTolerance',1e-18,'MaxFunEvals',10000,'StepTolerance', 1e-18, 'OptimalityTolerance', 1e-18,'MaxIterations', 10000); %make sure the solver doesn't end too early
    c1=lsqcurvefit(f1, c0, t, i, lb, ub, options);
    
    RMSE=sqrt(sum(((i-f1(c1,t)).^2))./length(i));
    
    if k==10
    
    figure
    subplot(2,1,1)
    plot((t-0.1),i,'-k',(t-0.1),f1(c1,t),'r-.')
    legend('Experimental Data','Fitted Data')
    no=xlabel('Time (s)');
    ye=ylabel('Current (A/g)');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    
    subplot(2,1,2)
    plot((t-0.1), (i-f1(c1,t)), '-k')
    grid on
    no=xlabel('Time (s)');
    ye=ylabel('Residuals');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_fitting10.fig');
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_fitting10.bmp');
   
    end
    
    if k==16
    
    figure
    subplot(2,1,1)
    plot((t-0.1),i,'-k',(t-0.1),f1(c1,t),'r-.')
    legend('Experimental Data','Fitted Data')
    no=xlabel('Time (s)');
    ye=ylabel('Current (A/g)');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    
    subplot(2,1,2)
    plot((t-0.1), (i-f1(c1,t)), '-k')
    grid on
    no=xlabel('Time (s)');
    ye=ylabel('Residuals');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_fitting16.fig');
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_fitting16.bmp');
    
    figure
    plot((t-0.1), (E/c1(1)).*exp(-t./(c1(1)*c1(2)))+(E/c1(3)).*exp(-t./(c1(3)*c1(4)))+c1(5)./t.^0.5+c1(6)*ones(length(t),1), '-k', (t-0.1), (E/c1(1)).*exp(-t./(c1(1)*c1(2))), '-.k', (t-0.1), (E/c1(3)).*exp(-t./(c1(3)*c1(4))), '--k', (t-0.1), c1(5)./t.^0.5,':k', (t-0.1), c1(6)*ones(length(t),1), '-.r')
    no=xlabel('Time (s)');
    ye=ylabel('Current (A/g)');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    legend('i_{t}', 'i_{DL1}','i_{DL2}', 'i_D', 'i_R')
    xlim([0 15])
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_ibreakdown.fig');
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_ibreakdown.bmp');
    
    end
    
    fwd(k,:)=[c1(1),c1(2),c1(3),c1(4),c1(5),c1(6), RMSE];

end

Rs1=fwd(:,1); %assign columns to variables
Cdl1=fwd(:,2);
Rs2=fwd(:,3);
Cdl2=fwd(:,4);
Bc=fwd(:,5);
F=fwd(:,6);
Errorc=fwd(:,7);

%Now let us construct voltammograms for iD (both cumulative capacitance and
%total capacitance are covered here

Rs1=smooth(Rs1,'rlowess');
Cdl1=smooth(Cdl1,'rlowess');
Rs2=smooth(Rs2,'rlowess');
Bc=smooth(Bc,'rlowess');
Cdl2=smooth(Cdl2,'rlowess');
F=smooth(F,'rlowess');
iD=zeros(fna,1);

Potential=linspace(v_min,v_max,fna);
Potentialr=linspace((v_max-0.025),(v_min+0.025),rn);
srate=(1000*E)./linspace(0.1,300,3000);
timeer=linspace(0.1,300,3000);
savd=zeros(3000,fna);

for j=1:fna % This section calculates stepwise voltammograms for iD
    iDa=Bc(j,1)./(linspace(0.1,300,3000).^0.5);
    iDb=Bc(j,1)./(linspace(0.1,time,time*5).^0.5);
    savd(:,j)=iDa';
    
    for p=1:3000
    iD(p,j)=(mean(savd(1:p,j)));
       
        if j<2
        Charged(p,j)=iD(p,j)*timeer(1,p);
        elseif j>1
        Charged(p,j)=((iD(p,j))*timeer(1,p))+Charged(p,j-1);
        end
        
    end
    
end

for j=1:fna % This section calculates stepwise voltammograms for iD
    iDex=Bc(j,1)./(linspace(0.2,time,time*5).^0.5);
    iDe(j,:)=mean(iDex);

end

temp=cumtrapz(Potential, iDe);
temp1=temp'-[0 temp(1:(fna-1))'];
temp1(1,1)=temp(2,1)-(temp(2,1)/20);
ap=temp1./(0.025^2/time);
savdl1=zeros(3000,fna);

for j=1:fna % This section calculates stepwise voltammograms for iD
    idl1=(E/Rs1(j,1)).*exp(-timeer./(Rs1(j,1)*Cdl1(j,1)));
    savdl1(:,j)=idl1';

for p=1:3000
    iDl1(p,j)=(mean(savdl1(1:p,j)));
        
        if j<2
        Chargedl1(p,j)=iDl1(p,j)*timeer(1,p);
        elseif j>1
        Chargedl1(p,j)=((iDl1(p,j))*timeer(1,p))+Chargedl1(p,j-1);
        end
        
end

end

savdl2=zeros(3000,fna);

for j=1:fna % This section calculates stepwise voltammograms for iD
    idl2=(E/Rs2(j,1)).*exp(-timeer./(Rs2(j,1)*Cdl2(j,1)));
    savdl2(:,j)=idl2';
    
    for p=1:3000
    iDl2(p,j)=mean(savdl2(1:p,1));
        
        if j<2
        Chargedl2(p,j)=iDl2(p,j)*timeer(1,p);
        elseif j>1
        Chargedl2(p,j)=((iDl2(p,j))*timeer(1,p))+Chargedl2(p,j-1);
        end
        
    end
    
end

Fful=ones(3000,1)*F';

for j=1:fna % This section calculates stepwise voltammograms for iD
   
    for p=1:3000
        
        if j<2
        Chargeff(p,j)=Fful(p,j)*timeer(1,p);
        elseif j>1
        Chargeff(p,j)=((Fful(p,j))*timeer(1,p))+Chargeff(p,j-1);
        end
        
    end
    
end

%Now for the reverse cycle
Er=E;
rvs=zeros(rn,7);

for z=1:rn                        
    D=C(3000*z+1:3000*z+time*5,:);
    t=(linspace(0.1, time, time*5))';
    i=-1*i2(3000*z+1:3000*z+time*5,:);
    lb=[0 0 0 0 0 -5]; % lower bound on constants (whatever is reasonable)
    ub=[10 400 10 200 20 5]; % upper bound on constants
    c0=[0.05 100 0.0020 20 1 0]; %Initial estimates (as close as possible to true values, use EIS results etc. if necessary, can get some very unlikely results otherwise)
    fn=@(c,t) (Er/c(1))*exp(-t./(c(1)*c(2)))+(Er/c(3))*exp(-t./(c(3)*c(4)))+(c(5)./(t.^0.5))+c(6);
    options = optimoptions('lsqcurvefit','FunctionTolerance',1e-18,'MaxFunEvals',10000,'StepTolerance', 1e-18, 'OptimalityTolerance', 1e-18,'MaxIterations', 10000); %make sure the solver doesn't end too early
    c1=lsqcurvefit(fn, c0, t, i, lb, ub, options);
    RMSE=sqrt(sum(((i-fn(c1,t)).^2))./length(i));
    
    if z==10
    figure
    subplot(2,1,1)
    plot((t-0.1),i,'-k',(t-0.1),fn(c1,t),'r-.')
    legend('Experimental Data','Fitted Data')
    xlabel('Time (s)')
    ylabel('Current (A/g)')
    no=xlabel('Time (s)');
    ye=ylabel('Current (A/g)');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    subplot(2,1,2)
    plot((t-0.1), (i-fn(c1,t)), '-k')
    grid on
    no=xlabel('Time (s)');
    ye=ylabel('Residuals');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_2_fitting10.fig');
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_2_fitting10.bmp');
    SSR10=1000*sum((i-fn(c1,t)).^2);
    
    end
    
    if z==30
    figure
    subplot(2,1,1)
    plot((t-0.1),i,'-k',(t-0.1),fn(c1,t),'r-.')
    legend('Experimental Data','Fitted Data')
    xlabel('Time (s)')
    ylabel('Current (A/g)')
    no=xlabel('Time (s)');
    ye=ylabel('Current (A/g)');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    subplot(2,1,2)
    plot((t-0.1), (i-fn(c1,t)), '-k')
    grid on
    no=xlabel('Time (s)');
    ye=ylabel('Residuals');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_2_fitting30.fig');
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_2_fitting30.bmp');
    figure
    plot((t-0.1), (Er/c1(1)).*exp(-t./(c1(1)*c1(2)))+(Er/c1(3)).*exp(-t./(c1(3)*c1(4)))+c1(5)./t.^0.5+c1(6)*ones(length(t),1), '-k', (t-0.1), (Er/c1(1)).*exp(-t./(c1(1)*c1(2))), '-.k', (t-0.1), (Er/c1(3)).*exp(-t./(c1(3)*c1(4))), '--k', (t-0.1), c1(5)./t.^0.5,':k', (t-0.1), c1(6)*ones(length(t),1), '-.r')
    no=xlabel('Time (s)');
    ye=ylabel('Current (A/g)');
    set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    legend('i_{t}', 'i_{DL1}','i_{DL2}', 'i_D', 'i_R')
    xlim([0 15])
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_ibreakdownrvs.fig');
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_ibreakdownrvs.bmp');
    
    end
    rvs(z,:)=[c1(1),c1(2),c1(3),c1(4),c1(5),c1(6), RMSE];

end

Rs1r=rvs(:,1); %assign columns to variables
Cdl1r=rvs(:,2);
Rs2r=rvs(:,3);
Cdl2r=rvs(:,4);
Br=-1*rvs(:,5);
Fr=-1*rvs(:,6);
Errorcr=rvs(:,7);

Rs1r=smooth(Rs1r,'rlowess');
Cdl1r=smooth(Cdl1r,'rlowess');
Rs2r=smooth(Rs2r,'rlowess');
Cdl2r=smooth(Cdl2r,'rlowess');
Fr=smooth(Fr,'rlowess');
Br=smooth(Br,'rlowess');
tan=linspace(300, max(T1), fna);
tca=linspace(300,max(T2),rn);
t_tr=[tan tca+max(T1)];
Rs1_tr=[Rs1.' Rs1r.'];
Rs2_tr=[Rs2.' Rs2r.'];
Cdl1_tr=[Cdl1.' Cdl1r.']; 
Cdl2_tr=[Cdl2.' Cdl2r.'];
CDL=Cdl1_tr+Cdl2_tr;
avCdl=mean(Cdl1_tr+Cdl2_tr);
avRs1=mean(Rs1_tr);
avRs2=mean(Rs2_tr);

for x=1:rn
    iDexr=Br(x,1)./((linspace(0.1,time,time*5)).^0.5);
    ire(x,:)=mean(iDexr);
end

tempr=cumtrapz(Potentialr, abs(ire));
tempr1=tempr'-[0 tempr(1:(rn-1))'];
tempr1(1,1)=tempr(2,1)-(tempr(2,1)/20);
apr=tempr1./(0.025^2/time);

%plot a sequential comparison of various R and C values

figure
hold on
subplot(4,1,1)
yyaxis left
a1=plot((t_tr-0.1), Rs1_tr, '-.k');
no=xlabel('Time (s)');
ye=ylabel('Resistance (\Omega\cdotg)');
set(gca, 'FontName', 'Arial', 'Fontsize', 10);
set(no,'Fontname', 'Arial', 'Fontsize', 10);
set(ye, 'Fontname', 'Arial', 'Fontsize', 10);
ylim([min(Rs1_tr)-0.005 max(Rs1_tr)+0.005])
set(gca, 'YColor', 'k') %set colour of axis to black
text(0.02,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
box on
text(10000, max(Rs1_tr)+0.01, 'Potential (V vs SCE)','HorizontalAlignment','Center','Fontsize', 10)
text(0,max(Rs1_tr)+0.007,'-0.5','HorizontalAlignment','Center','Fontsize', 10)
text(5000, max(Rs1_tr)+0.007,'-0.1','HorizontalAlignment','Center','Fontsize', 10);
text(10000,max(Rs1_tr)+0.007,'0.3','HorizontalAlignment','Center','Fontsize', 10);
text(15000,max(Rs1_tr)+0.007,'-0.1','HorizontalAlignment','Center','Fontsize', 10);
text(20000,max(Rs1_tr)+0.007,'-0.5','HorizontalAlignment','Center','Fontsize', 10);
ax1 = gca; % current axes
SP=10000; %your point goes here 
line([SP SP],get(ax1,'YLim'),'Color',[1 0 0], 'Linestyle', ':')

yyaxis right
a2=plot(t_tr-0.1,Cdl1_tr,'-k');
ze=ylabel('Capacitance (F/g)');
set(ze, 'Fontname', 'Arial', 'Fontsize', 10)
set(gca, 'YColor', 'k')
legend([a1,a2],{'R_{S1}', 'C_{DL1}'},'Location', 'eastoutside');

subplot(4,1,2)
yyaxis left
a1=plot((t_tr-0.1), Rs2_tr, '-.k');
no=xlabel('Time (s)');
ye=ylabel('Resistance (\Omega\cdotg)');
set(gca, 'FontName', 'Arial', 'Fontsize', 10);
set(no,'Fontname', 'Arial', 'Fontsize', 10);
set(ye, 'Fontname', 'Arial', 'Fontsize', 10);
ylim([min(Rs2_tr)-0.002 max(Rs2_tr)+0.002])
set(gca, 'YColor', 'k') %set colour of axis to black
text(0.02,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
box on
text(10000, max(Rs2_tr)+0.004, 'Potential (V vs SCE)','HorizontalAlignment','Center','Fontsize', 10)
text(0,max(Rs2_tr)+0.003,'-0.5','HorizontalAlignment','Center','Fontsize', 10)
text(5000, max(Rs2_tr)+0.003,'-0.1','HorizontalAlignment','Center','Fontsize', 10);
text(10000,max(Rs2_tr)+0.003,'0.3','HorizontalAlignment','Center','Fontsize', 10);
text(15000,max(Rs2_tr)+0.003,'-0.1','HorizontalAlignment','Center','Fontsize', 10);
text(20000,max(Rs2_tr)+0.003,'-0.5','HorizontalAlignment','Center','Fontsize', 10);
ax1 = gca; % current axes
SP=10000; %your point goes here 
line([SP SP],get(ax1,'YLim'),'Color',[1 0 0], 'Linestyle', ':')

yyaxis right
a2=plot(t_tr-0.1,Cdl2_tr,'-k');
ze=ylabel('Capacitance (F/g)');
set(ze, 'Fontname', 'Arial', 'Fontsize', 10)
set(gca, 'YColor', 'k')
legend([a1,a2],{'R_{S2}', 'C_{DL2}'},'Location', 'eastoutside');

subplot(4,1,3)
frac1=[Cdl1' Cdl1r']./[Cdl1'+Cdl2'+ap Cdl1r'+Cdl2r'+(-1*apr)];
frac2=[Cdl2' Cdl2r']./[Cdl1'+Cdl2'+ap Cdl1r'+Cdl2r'+(-1*apr)];
frac3=[ap (-1*apr)]./[Cdl1'+Cdl2'+ap Cdl1r'+Cdl2r'+(-1*apr)];
plot((t_tr-0.1), 100*frac1,'-.k',(t_tr-0.1),100*frac2,'--k',(t_tr-0.1),100*frac3,'-k')
no=xlabel('Time (s)');
ye=ylabel('Fractional Contribution (%)');
set(gca, 'FontName', 'Arial', 'Fontsize', 10);
set(no,'Fontname', 'Arial', 'Fontsize', 10);
set(ye, 'Fontname', 'Arial', 'Fontsize', 10);
set(gca, 'YColor', 'k') %set colour of axis to black
ylim([0 100])
text(10000, 115, 'Potential (V vs SCE)','HorizontalAlignment','Center','Fontsize', 10)
text(0,105,'-0.5','HorizontalAlignment','Center','Fontsize', 10)
text(5000, 105,'-0.1','HorizontalAlignment','Center','Fontsize', 10);
text(10000,105,'0.3','HorizontalAlignment','Center','Fontsize', 10);
text(15000,105,'-0.1','HorizontalAlignment','Center','Fontsize', 10);
text(20000,105,'-0.5','HorizontalAlignment','Center','Fontsize', 10);
ax1 = gca; % current axes
SP=10000; %your point goes here 
line([SP SP],get(ax1,'YLim'),'Color',[1 0 0], 'Linestyle', ':')
text(0.02,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
box on
legend('C_{DL1}','C_{DL2}','C_D','Location', 'eastoutside')

subplot(4,1,4)
hold on
%start with the left axis plot (requires matlab 2016 edition or later for yyaxis command to work)
plot((t_tr-0.1), [Cdl1'+Cdl2'+ap Cdl1r'+Cdl2r'+(-1*apr)],'-k',(t_tr-0.1), [Cdl1' Cdl1r'],'-.k', (t_tr-0.1),[Cdl2' Cdl2r'],':k',(t_tr-0.1), [ap (-1*apr)], '--k')
no=xlabel('Time (s)');
ye=ylabel('Capacitance (F/g)');
set(gca, 'FontName', 'Arial', 'Fontsize', 10);
set(no,'Fontname', 'Arial', 'Fontsize', 10);
set(ye, 'Fontname', 'Arial', 'Fontsize', 10)
set(gca, 'YScale', 'log')

set(gca, 'YColor', 'k') %set colour of axis to black
text(0.02,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
box on
text(10000, 1500, 'Potential (V vs SCE)','HorizontalAlignment','Center','Fontsize', 10,'Fontname','Arial')
text(0,1250,'-0.5','HorizontalAlignment','Center','Fontsize', 10,'Fontname','Arial')
text(5000, 1250,'-0.1','HorizontalAlignment','Center','Fontsize', 10,'Fontname','Arial');
text(10000,1250,'0.3','HorizontalAlignment','Center','Fontsize', 10,'Fontname','Arial');
text(15000,1250,'-0.1','HorizontalAlignment','Center','Fontsize', 10,'Fontname','Arial');
text(20000,1250,'-0.5','HorizontalAlignment','Center','Fontsize', 10,'Fontname','Arial');
plot([10000 10000], [0.1 1000],':r')
hlegend=legend('C_{Total}', 'C_{DL1}', 'C_{DL2}','C_D', 'Location', 'best');
hlegend.Location='eastoutside';
xlim([0 19500]);
ylim([10 1000]);
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_subplots.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_subplots.bmp');

clear ir
ir=zeros(rn,1);
savdr=zeros(3000,rn);

for j=1:rn
    iDr=Br(j,1)./((linspace(0.1,300,3000)).^0.5);
    savdr(:,j)=iDr';
    
    for p=1:3000
    ir(p,j)=mean(savdr(1:p,j));
        
        if j<2
        Chargedr(p,j)=ir(p,j)*timeer(1,p);
        elseif j>1
        Chargedr(p,j)=((ir(p,j))*timeer(1,p))+Chargedr(p,j-1);
        end
        
    end
    
end

savdl1r=zeros(3000,rn);

for j=1:rn % This section calculates stepwise voltammograms for iD
    idl1r=(E/Rs1r(j,1)).*exp(-timeer./(Rs1r(j,1)*Cdl1r(j,1)));
    savdl1r(:,j)=idl1r';
    
    for p=1:3000
    iDl1r(p,j)=mean(savdl1r(1:p,j));
        
        if j<2
        Chargedl1r(p,j)=iDl1r(p,j)*timeer(1,p);
        elseif j>1
        Chargedl1r(p,j)=((iDl1r(p,j))*timeer(1,p))+Chargedl1r(p,j-1);
        end
        
    end
    
end

savdl2r=zeros(3000,rn);

for j=1:rn % This section calculates stepwise voltammograms for iD
    idl2r=(E/Rs2r(j,1)).*exp(-timeer./(Rs2r(j,1)*Cdl2r(j,1)));
    savdl2r(:,j)=idl2r';
    
    for p=1:3000
    iDl2r(p,j)=mean(savdl2r(1:p,j));
        
        if j<2
        Chargedl2r(p,j)=iDl2r(p,j)*timeer(1,p);
        elseif j>1
        Chargedl2r(p,j)=((iDl2r(p,j))*timeer(1,p))+Chargedl2r(p,j-1);
        end
        
    end
    
end

Ffulr=ones(3000,1)*Fr';

for j=1:rn % This section calculates stepwise voltammograms for iD
    
    for p=1:3000
        
        if j<2
        Chargeffr(p,j)=Ffulr(p,j)*timeer(1,p);
        elseif j>1
        Chargeffr(p,j)=((Ffulr(p,j))*timeer(1,p))+Chargeffr(p,j-1);
        end
        
    end
    
end

for p=1:3000
    Anocapd(1,p)=(max([Charged(p,1:end) Chargedr(p,1:end)])-min(Charged(p,1:end)))/0.8;
end

for p=1:3000
    Anocapdl1(1,p)=(max([Chargedl1(p,1:end) Chargedl1r(p,1:end)])-min(Chargedl1(p,1:end)))/0.8;
end

for p=1:3000
    Anocapdl2(1,p)=(max([Chargedl2(p,1:end) Chargedl2r(p,1:end)])-min(Chargedl2(p,1:end)))/0.8;
end

for p=1:3000
    Anocapff(1,p)=(max([Chargeff(p,1:end) Chargeffr(p,1:end)])-min(Chargeff(p,1:end)))/0.8;
end

SEad=((Anocapd.*(0.8^2)./2)/3.6);
SEadl1=((Anocapdl1.*(0.8^2)./2)/3.6);
SEadl2=((Anocapdl2.*(0.8^2)./2)/3.6);
SEaff=((Anocapff.*(0.8^2)./2)/3.6);
SEc=SEad+SEadl1+SEadl2;
SPPad=SEad./((0.8./(srate./1000))./3600);
SPPadl1=SEadl1./((0.8./(srate./1000))./3600);
SPPadl2=SEadl2./((0.8./(srate./1000))./3600);
SPPaff=SEaff./((0.8./(srate./1000))./3600);
SPPc=SPPad+SPPadl1+SPPadl2;
f1x=[SEc(end) SEad(end) SEadl1(end) SEadl2(end)];
f1y=[SPPc(end) SPPad(end) SPPadl1(end) SPPadl2(end)];
f2x=[SEc(1,125) SEad(1,125) SEadl1(1,125) SEadl2(1,125)];
f2y=[SPPc(1,125) SPPad(1,125) SPPadl1(1,125) SPPadl2(1,125)];
f3x=[SEc(1,10) SEad(1,10) SEadl1(1,10) SEadl2(1,10)];
f3y=[SPPc(1,10) SPPad(1,10) SPPadl1(1,10) SPPadl2(1,10)];
f4x=[SEc(1,1) SEad(1,1) SEadl1(1,1) SEadl2(1,1)];
f4y=[SPPc(1,1) SPPad(1,1) SPPadl1(1,1) SPPadl2(1,1)];

for p=1:3000
    Ccapd(1,p)=(max([Charged(p,1:end) Chargedr(p,1:end)])-Charged(p,end)+min(Charged(p,1:end)))/0.8;
end

for p=1:3000
    Ccapdl1(1,p)=(max([Chargedl1(p,1:end) Chargedl1r(p,1:end)])-Chargedl1(p,end)+min(Chargedl1(p,1:end)))/0.8;
end

for p=1:3000
   Ccapdl2(1,p)=(max([Chargedl2(p,1:end) Chargedl2r(p,1:end)])-Chargedl2(p,end)+min(Chargedl2(p,1:end)))/0.8;
end

for p=1:3000
    Ccapff(1,p)=(max([Chargeff(p,1:end) Chargeffr(p,1:end)])-Chargeff(p,end)+min(Chargeff(p,1:end)))/0.8;
end

SEcd=((Ccapd.*(0.8^2)./2)/3.6);
SEcdl1=((Ccapdl1.*(0.8^2)./2)/3.6);
SEcdl2=((Ccapdl2.*(0.8^2)./2)/3.6);
SEcff=((Ccapff.*(0.8^2)./2)/3.6);
SPPcd=SEcd./((0.8./(srate./1000))./3600);
SPPcdl1=SEcdl1./((0.8./(srate./1000))./3600);
SPPcdl2=SEcdl2./((0.8./(srate./1000))./3600);
SPPcff=SEcff./((0.8./(srate./1000))./3600);

Fcap=F.*300/0.8;
Frcap=Fr.*300/0.8;
Bcap=(2.*Bc.*sqrt(300))./0.8;
Brcap=(2.*Br.*sqrt(300))./0.8;
Cdl1cap=(E.*Cdl1.*exp(-300./Cdl1.*Rs1)-1)./0.8;
Cdl1rcap=(-E.*Cdl1r.*exp(-300./Cdl1r.*Rs1r)-1)./0.8;
Cdl2cap=(E.*Cdl2.*exp(-300./Cdl2.*Rs2)-1)./0.8;
Cdl2rcap=(-E.*Cdl2r.*exp(-300./Cdl2r.*Rs2r)-1)./0.8;

%Begin compiling variables

FD=[F' Fr' F(1,1)];
Potentiala=[Potential Potentialr Potential(1,1)];
Rs1_tr=[Rs1.' Rs1r.'];
Rs2_tr=[Rs2.' Rs2r.'];
Cdl1_tr=[Cdl1.' Cdl1r.']; 
Cdl2_tr=[Cdl2.' Cdl2r.'];
CDL=Cdl1_tr+Cdl2_tr;
avCdl=mean(Cdl1_tr+Cdl2_tr);
avRs1=mean(Rs1_tr);
avRs2=mean(Rs2_tr);
clear x
clear y
x=linspace(0,19500,195000);
y=-0.525+0.025*(heaviside(x-300)+heaviside(x-600)+heaviside(x-900)+heaviside(x-1200)+heaviside(x-1500)+heaviside(x-1800)+heaviside(x-2100)+heaviside(x-2400)+heaviside(x-2700)+heaviside(x-3000)+heaviside(x-3300)+heaviside(x-3600)+heaviside(x-3900)+heaviside(x-4200)+heaviside(x-4500)+heaviside(x-4800)+heaviside(x-5100)+heaviside(x-5400)+heaviside(x-5700)+heaviside(x-6000)+heaviside(x-6300)+heaviside(x-6600)+heaviside(x-6900)+heaviside(x-7200)+heaviside(x-7500)+heaviside(x-7800)+heaviside(x-8100)+heaviside(x-8400)+heaviside(x-8700)+heaviside(x-9000)+heaviside(x-9300)+heaviside(x-9600)+heaviside(x-9900)-(heaviside(x-10200)+heaviside(x-10500)+heaviside(x-10800)+heaviside(x-11100)+heaviside(x-11400)+heaviside(x-11700)+heaviside(x-12000)+heaviside(x-12300)+heaviside(x-12600)+heaviside(x-12900)+heaviside(x-13200)+heaviside(x-13500)+heaviside(x-13800)+heaviside(x-14100)+heaviside(x-14400)+heaviside(x-14700)+heaviside(x-15000)+heaviside(x-15300)+heaviside(x-15600)+heaviside(x-15900)+heaviside(x-16200)+heaviside(x-16500)+heaviside(x-16800)+heaviside(x-17100)+heaviside(x-17400)+heaviside(x-17700)+heaviside(x-18000)+heaviside(x-18300)+heaviside(x-18600)+heaviside(x-18900)+heaviside(x-19200)+heaviside(x-19500)+heaviside(x-19800)+heaviside(x-20400)));

%___________________________________________________________________________________
%Here we make our major plots

figure
hold on
yyaxis left %start with the left axis plot (requires matlab 2016 edition or later for yyaxis command to work)
plot((t_tr-0.1), Rs1_tr, '-k.', (t_tr-0.1), Rs2_tr, '-.k');
set(gca,'YScale','log')
xlim([0 19500]);
ylim([0 0.1]);
no=xlabel('Time (s)');
ye=ylabel('Resistance (\Omega\cdotg)');
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
set(no,'Fontname', 'Arial', 'Fontsize', 14);
set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
set(gca, 'YColor', 'k') %set colour of axis to black
box on

yyaxis right %onwards to the right axis plot (Potential)
plot(x, y, 'r');
ze=ylabel('Potential (V vs SCE)');
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
set(ze, 'Fontname', 'Arial', 'Fontsize', 14)
set(gca, 'YColor', 'k')
alegend=legend('R_{S1}', 'R_{S2}','Potential', 'Location', 'northwest');
hold off
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_Rvt.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_Rvt.bmp');

%Plot the cathodic capacitances
figure
semilogx(srate, Anocapd+Anocapdl1+Anocapdl2, '-k',...
    srate,Anocapdl1, '-.k'...
    ,srate,Anocapdl2, '--k',...
    srate,Anocapd,':k')
legend('C_{Total}','C_{DL1}', 'C_{DL2}','C_{D}');
no=xlabel('Scan Rate (mV/s)');
ye=ylabel('Specific Capacitance (F/g)');
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
set(no,'Fontname', 'Arial', 'Fontsize', 14);
set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
xlim([0.003 125])
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_SRcapacitance.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_SRCapacitiance.bmp');

figure
loglog(SEad+SEadl1+SEadl2, SPPad+SPPadl1+SPPadl2, '-k',SEadl1 , SPPadl1, '--k', SEadl2, SPPadl2, ':k', SEad, SPPad, '-.k',f1x,f1y,':r',f2x,f2y,':r',f3x,f3y,':r',f4x,f4y,':r')
no=xlabel('Specific Energy (Wh/kg)');
ye=ylabel('Specific Power (W/kg)');
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
set(no,'Fontname', 'Arial', 'Fontsize', 14);
set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
lgndc=legend('Overall Performance', 'C_{DL1} Effects', 'C_{DL2} Effects', 'C_{D} Effects');
lgndc.Location='southeast';
text(SEc(end), SPPc(end)+10, '0.083 mV/s','FontName','Arial', 'Fontsize',9)
text(SEc(1,125),SPPc(1,125)+15,'0.1 mV/s','FontName','Arial', 'Fontsize',9)
text(SEc(1,10),SPPc(1,10)+100, '12.5 mV/s','FontName','Arial', 'Fontsize',9)
text(SEc(1,1),SPPc(1,1)+700,'125 mV/s','FontName','Arial', 'Fontsize',9)
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_Rag.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_Rag.bmp');

figure %Plot the cathodic ragone diagram alone
loglog(SEad+SEadl1+SEadl2, SPPad+SPPadl1+SPPadl2, '-k')
no=xlabel('Specific Energy (Wh/kg)');
ye=ylabel('Specific Power (W/kg)');
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
set(no,'Fontname', 'Arial', 'Fontsize', 14);
set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
lgndc=legend('Overall Performance');
lgndc.Location='Best';
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_Ragalone.fig');
saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\SPECS Results\Copper\Cu_1_org_2_1_Ragalone.bmp');











