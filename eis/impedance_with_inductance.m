%____________________________________
%
% Author        - Vitamin-C
%
% Status        - Functional
%
% Description   - Models and extracts multiple capacitance, inductive, and
%                 resistance values across a set potential range for supercapacitve
%                 materials. Can be modified for battery-like behaviour.
% 
% Use Comments  - Adjust the newby file (which contains the mathematical
%                 model of the impedence spectra of our material, which
%                 exhibits inductance, double-layer, and diffusion effects.
%                 Set up save dir for first time use
%                 Use find and replace to change variable fields.
%                 This script could be improved with descriptive variable
%                 names.

clear all
clc

% Dependancies - newby.m

% Set up variables and parameters
% File - Ti_4
% Save Folder - EIS\Untreated
v_min=-0.5; %lower limit on potential
v_max=0.3; %upper limit on potential
n_s=((v_max-v_min)/0.025)+1; % number of scans, dependant on potential window
n_s=round(n_s);

%_______________________________________________________

data=zeros(n_s,14); %this is an array of zeros in which to slot our data, the first term is the rows, second is columns, in this case EIS was taken in 33 files, I have 8 variables and wish to record Capacitance as well
Potential=linspace(v_min, v_max, n_s);

 %______________________________________________________________________________
 
 for k=1:n_s % loop to create temporary names for each file
    myfilename=sprintf('%s_%d','C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\Capacitance stuff\Ti_4\EIS_1\Ti_4_',k); %set filenames, directory goes here
   
    %____________________________________________________________________________________________________________________
    m=0.00089196; %scaled sample mass here (usually SSA is used, but mass is simpler and I prefer F/g to F/cm2)

    %____________________________________________________________________________________________________________________
    A=importdata(myfilename, ' ');%Import as a matrix delimited by the second term (a space for this data)
    A=A.data; %Get rid of those pesky heading lines
    time=A(:,1); %separate and adjust variables
    Zr=A(:,2).*m; %modelling is problematic for small values, so we scale here to milliOhm/g
    Zi=A(:,3).*m;
    Zim=(-1*Zi);
    Hz=A(:,4);
    rads=2*pi*Hz;
    Zcom=sqrt(((Zr.*m).^2)+((Zi.*m).^2)); %This is the magnitude of the complex electric impedance (in Ohm.g)
    Cf=abs((1./(rads.*Zim))); %This is the frequency dependant capacitance of the material (F/g)
    Cfmax=Cf(28,1); %The highest frequency dependant capacitance outside the impedance region
   %____________________________________________________________________________________________________________
    
   for j=1:20                  % Tells us where the circle turns and ends, so is used to estimate initial parameters
       
       if Zim(j,1)>Zim(j+1,1)
           wh(1,j)=Zr(j,1);
           wa(1,j)=Zim(j,1);
           fa(1,j)=Hz(j,1);
       end
       
   end
   
   wht=wh';                           % this part gets rid of zeros from the for loop
   idx=wht==0;
   bwh(~sort(idx))=wht(~idx);   
    %___________________________________________________________________________________________
    
    Rsest=Zr(1,1);
    Rctest=0.8;
    s1est=25000;
    s2est=0.1;
    m1est=0.9;                     %Estimates of variables
    m2est=0.9; %how 'circular' the circle is
    fc=max(fa);  %This is the characteristic frequency for a potential, should be the same for all frequencies
    
    %__________________________________________________________________
    
   c0=[6.0379e-03    4.0885e-02    8.2451e-02    1.7186e-02    -3.0729e-02    1.5696e-01    -6.7493e-01    2.3305e-02    7.2047e-01];  %This is the model itself, it relies on files f1 to f13 which contain functions, need to have them either open or in the same directory
    lb=[0,0,0,0,-1000,0,-1,0,0];
    ub=[10,50,500000,1,0,20,0,100,1];
    options = optimoptions('fmincon', 'MaxFunctionEvaluations',100000, 'MaxIterations', 100000, 'ConstraintTolerance', 1e-20, 'Algorithm', 'interior-point','StepTolerance',1e-20, 'FunctionTolerance', 1e-20); 
    format short e % output c vector in scientific notation
    f=@(c) newby(c, rads, Zr, Zim);
    c1=fmincon(f, c0,[],[],[],[],lb, ub,[],options); %For good data this will be quick and accurate, for less neat data not so much
    
    %__________________________________________________________________PLOTS
    
    Cap=((((c1(2))/(c1(3)))^(1/c1(4)))/(c1(2))); %an alternate method of calculating capacitance, it is low as it is only on a small potential perturbation (about 10 mV)
    f_k=figure;box on;
    subplot(2,1,1)
    plot(Zr, Zim, '-k', Zrt(c1,rads), -Zit(c1,rads), ':r') %This will show if the model has produced a good fit
    no=xlabel('Z_{Re} (\Omega\cdotg)');
    ye=ylabel('-Z_{Im} (\Omega\cdotg)');
    legend('Experimental Data', 'Fitted Data', 'Location','northwest')
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    xlim([0 max(Zr)])
    ylim([0 max(Zim)])
    
    subplot(2,1,2) %This subplot plots frequency - dependant capacitance against frequency
    plot(Hz,Cf,'-k')
    legend('Capacitance')
    no=xlabel('\nu (Hz)');
    ye=ylabel('C_\nu (F/g)');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    ax = gca;
    c = ax.XScale;
    ax.XScale = 'log';
    ax = gca;
    c = ax.YScale;
    ax.YScale = 'log';
    xlim([0.1 100000])
    ylim([0.1 1000])
    
    %___________________________________________________________________________________________________
    %Create a table of the variables and data
    data(k,:)=[Potential(1,k), c1(1),c1(2),c1(3),c1(4),c1(5),c1(6),c1(7),c1(8),c1(9),Cap,fc,Cfmax,f(c1)]; % Put the data in the array we set up at the top, this should have as many entries as you want to store
    
    % Save a couple of the plots so they can be used later
    
    if k==1
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(-0.5).fig');
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(-0.5).bmp');
    end
    
    if k==10
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(-0.275).fig');
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(-0.275).bmp');
    end
    
    if k==20
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(-0.025).fig');
     saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(-0.025).bmp');
    end
    
    if k==33
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(0.3).fig'); 
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Ti_4_1(-0.3).bmp');
    end
    
 end
 %__________________________________________________________________________________________________
 
 Pot=data(:,1); %Matlab does not allow subscript etc in variable names so need to change the tabulated names later
    R_s=data(:,2);
    Rct=data(:,3);
    s1=data(:,4);
    m1=data(:,5);
    Rl=data(:,6);
    L=data(:,7);
    m2=data(:,8);
    s3=data(:,9);
    m3=data(:,10);
    C=data(:,11);
    FC=data(:,12);
    CFm=data(:,13);
    SSE=(data(:,14)./100000);
    var=table(Pot,R_s,Rct,s1,m1,Rl,L,m2,s3,m3,C,FC,CFm,SSE); %Create and export a table of derived variable values
    my_directory='C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Tabulated data negative m2'; %Where to save goes here
    saveit  = fullfile( my_directory, 'Ti_4_1.csv' ); %Set a filename here
    writetable(var,saveit);
    
    %______________________________________________
    %Extra plots
    
    corr=C;
    corr(isnan(corr))=C(3,1);
    corr(1,1)=C(4,1);
    
    for k=2:n_s-3
       
       if C(k,1)>=5*C(k-1,1)
           corr(k,1)=((C(k-1,1)+C(k+1,1))/2);        
       end
       
    end
    
    pft=repnan(corr);
    pft=repnan(corr);
    linft=fit(Pot,pft,'Poly1');
    figure;box on
    hold on
    plot(Pot, C, '.k')
    plot(linft, '-r')
    no=xlabel('Potential (V vs SCE)');
    ye=ylabel('Capacitance (F/g)');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    legend('Capacitance', 'Location','northwest')
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsC_Ti_4_1.fig') %name for this one too
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsC_Ti_4_1.bmp');
    
    figure;box on
    plot(Pot,SSE, '.k')
    no=xlabel('Potential (V vs SCE)');
    ye=ylabel('S');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Errors_Ti_4_1.fig') %name for error plot
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\Errors_Ti_4_1.bmp') %name for error plot
  
    figure;box on
    plot(Pot, R_s, '.k')
    no=xlabel('Potential (V vs SCE)');
   ye=ylabel('Resistance (\Omega\cdotg)');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    legend('R_S', 'Location','best')
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsRs_Ti_4_1.fig') %name for this one too
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsRs_Ti_4_1.bmp') %name for this one too
    
 figure;box on
    plot(Pot, Rct, '.k')
    no=xlabel('Potential (V vs SCE)');
    ye=ylabel('Resistance (\Omega\cdotg)');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    legend('R_{CT}', 'Location','best')
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsRct_Ti_4_1.fig') %name for this one too
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsRct_Ti_4_1.bmp') %name for this one too
    
    
 figure;box on
 plot(Pot,L,'.k')
 no=xlabel('Potential (V vs SCE)');
    ye=ylabel('Inductance (H/g)');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    legend('L', 'Location','best')
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsL_Ti_4_1.fig') %name for this one too
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsL_Ti_4_1.bmp') %name for this one too
    
     figure;box on
 plot(Pot,Rl,'.k')
 no=xlabel('Potential (V vs SCE)');
    ye=ylabel('Resistance (\Omega\cdotg)');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    legend('R_{L}', 'Location','best')
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsRl_Ti_4_1.fig') %name for this one too
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsRl_Ti_4_1.bmp') %name for this one too
    
     figure;box on
    plot(Pot,CFm, '.k')
    no=xlabel('Potential (V vs SCE)');
    ye=ylabel('C_S^{max} (F/g)');
     set(gca, 'FontName', 'Arial', 'Fontsize', 12);
    set(no,'Fontname', 'Arial', 'Fontsize', 14);
    set(ye, 'Fontname', 'Arial', 'Fontsize', 14)
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsCsmax_Ti_4_1.fig') %name for error plot
    saveas(gcf, 'C:\Users\Fostiple\OneDrive - The University Of Newcastle\Documents\PhD\Experimental Records\Electrochemical Testing of Zee Carbons\Supercapacitor\EIS\Untreated\Modelled plots negative m2\PvsCsmax_Ti_4_1.bmp') %name for error plot
  