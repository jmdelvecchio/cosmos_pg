clear all
for m=1:5 %This is the number of MC iterations

clearvars -except m MC %This is to start each MC iteration fresh

%% Step 0: Choose your model domain parameters
Erate = 5.5; %erosion rate in m/Myr

rho_sap = 2.5; % rock density
rho_sed = 2.0; % fill density

pg_thickness = [5 3]'; %in m
periglacial_erosion = [0.5 1]'; %in m; syndepositional erosion
pg_deposit_times = [575 30]'; %in kya

%flag 0: base values
%flag 1: base values plus minus a normal distribution
%flag 2: a random number in some range

inherited_ratio_flag=1;
inherited_ratio_base=6.5; %applies in flags 0 and 1
inherited_ratio_range=0.5; %applies in flag 1
inherited_ratio_random_lower=5; %applies in flag 2
inherited_ratio_random_upper=9; %applies in flag 2

N10_inherited_flag=1;
N10_inherited_base=500000; %applies in flags 0 and 1
N10_inherited_range=100000; %applies in flag 1
N10_inherited_random_lower=0; %applies in flag 2
N10_inherited_random_upper=1000000; %applies in flag 2

%No uncertainty in production = 0;
%10% uncertainty in production = 1;
production_uncertainty_flag=1;

%%
pg_timestep = abs(diff([pg_deposit_times;[0]]))'; %just transforms
                                                  %user-entered times
Erate_cmyr = Erate.*(1e-4); %erosion rate in cm/yr
%% Step 1: numerically solve for steady-state profile

dz = 1; % depth increment in centimeters
zmax = 10000; % max depth is 4 m
z = (0:dz:zmax)';
MassDepth = z.*rho_sap;

dt = dz./Erate_cmyr; %time increment in years
t_max = 1000.*dt; %run time in years
t = (0:dt:t_max)';

if production_uncertainty_flag==1
    [P_total_10Be, P_total_26Al] = BMC_production_model_unc(MassDepth-0.5.*rho_sap.*dz); %simple exponential production for analytic solution - replace with something better
else
    [P_total_10Be, P_total_26Al] = BMC_production_model(MassDepth-0.5.*rho_sap.*dz); %simple exponential production for analytic solution - replace with something better
end

N_total_10Be = zeros(length(z),length(t)+1);
N_total_26Al = zeros(length(z),length(t)+1);

Be_halflife = 1387000;      % 10Be half life in years (Chmeleff et al., 2010)
Al_halflife = 705000;       % 26Al half life in years (Nishiizumi, 2004)
Be_decay_constant = log(2)/Be_halflife;
Al_decay_constant = log(2)/Al_halflife;


for i=1:length(t)
    for j = 1:(length(z)-1)
        N_total_10Be(j,i+1) = N_total_10Be(j+1,i) + P_total_10Be(j).*rho_sap.*dz.*dt...
            -Be_decay_constant.*N_total_10Be(j+1,i).*dt;
        N_total_26Al(j,i+1) = N_total_26Al(j+1,i) + P_total_26Al(j).*rho_sap.*dz.*dt...
            -Al_decay_constant.*N_total_26Al(j+1,i).*dt;         
    end
    
    N_total_10Be(length(z),i+1) = 0;
    N_total_26Al(length(z),i+1) = 0;

end

%% Step 2: erode top and copy only sample depth domain

N_total_10Be_paleosol=N_total_10Be(1:(19*100),end);
N_total_26Al_paleosol=N_total_26Al(1:(19*100),end);
rho_array = rho_sap.*ones(size(N_total_10Be_paleosol));

N_total_10Be_old=N_total_10Be_paleosol;
N_total_26Al_old=N_total_26Al_paleosol;

%% Step 3: add one or more packages of colluvium to the column 
for p=1:length(pg_thickness)
%% Step 3.1: erode any material syndepositionally
if periglacial_erosion(p) > 0 %New
    N_total_10Be_old = N_total_10Be_old((periglacial_erosion(p)*100)/dz:end);
    N_total_26Al_old = N_total_26Al_old((periglacial_erosion(p)*100)/dz:end);
    rho_array = rho_array((periglacial_erosion(p)*100)/dz:end);
end

%% Step 3.2 Create colluvial package
    %3.2.1 set inheritance ratio
    if inherited_ratio_flag == 0
       inherited_ratio = inherited_ratio_base; 
    elseif inherited_ratio_flag == 1
       inherited_ratio = inherited_ratio_base + inherited_ratio_range * randn(1,1);
    else
       inherited_ratio = inherited_ratio_random_lower + rand(1,1) * ...
           (inherited_ratio_random_upper - inherited_ratio_random_lower);
    end
    
    %3.2.2 N10 ratio
    if N10_inherited_flag == 0
        N10_inherited = N10_inherited_base * rho_sed * dz;
    elseif N10_inherited_flag == 1
        N10_inherited = (N10_inherited_base + N10_inherited_range*randn(1,1))...
            * rho_sed * dz;
    else
        N10_inherited = N10_inherited_random_lower +rand(1,1) * ...
            (N10_inherited_random_upper - N10_inherited_random_lower);
    end
    
    N26_inherited = N10_inherited*inherited_ratio; % inherited 26Al concentration (backcalculated)
    
    %% 3.2.3 deposit it
        N_total_10Be_pg = [(N10_inherited .* ones((pg_thickness(p)*100)/dz,1));N_total_10Be_old];
        N_total_26Al_pg = [(N26_inherited .* ones((pg_thickness(p)*100)/dz,1));N_total_26Al_old];
        rho_array = [(rho_sed .* ones((pg_thickness(p)*100)/dz,1));rho_array];
    z_temp = (0.1:dz:length(N_total_10Be_pg));

%% Step 4 Produce and decay nuclides for the given timelength
    for j=1:length(z_temp)
    N_total_10Be_cosmo(j) = N_total_10Be_pg(j) + P_total_10Be(j).*rho_array(j).*dz.*pg_timestep(p)*1000 ...
            -Be_decay_constant.*N_total_10Be_pg(j).*pg_timestep(p)*1000;
    
    N_total_26Al_cosmo(j) = N_total_26Al_pg(j) + P_total_26Al(j).*rho_array(j).*dz.*pg_timestep(p)*1000 ...
        -Al_decay_constant.*N_total_26Al_pg(j).*pg_timestep(p)*1000;         
    end
    N_total_10Be_old = N_total_10Be_cosmo';
    N_total_26Al_old = N_total_26Al_cosmo';
    
end %periglacial event loops

%% Step 5: Store final calculations in MC array at end of iteration
N_total_10Be = N_total_10Be_old;
N_total_26Al = N_total_26Al_old;
z = z_temp';

C_10Be = N_total_10Be./rho_array./dz;
C_26Al = N_total_26Al./rho_array./dz;
Ratio_26Al_10Be = C_26Al./C_10Be;

C_10Be_end = C_10Be(:,end);
C_26Al_end = C_26Al(:,end);

MC.C_10Be_mc(:,m)=C_10Be_end;
MC.C_26Al_mc(:,m)=C_26Al_end;
MC.Ratio_mc(:,m)=C_26Al_end./C_10Be_end;
end %MC loop

%% Step 6 calculate means and standard deviations of MC iterations
C_10Be_m = mean(MC.C_10Be_mc,2);
C_10Be_std = std(MC.C_10Be_mc,0,2);

C_26Al_m = mean(MC.C_26Al_mc,2);
C_26Al_std = std(MC.C_26Al_mc,0,2);

Ratio_m = mean(MC.Ratio_mc,2);
Ratio_std = std(MC.Ratio_mc,0,2);
%% Step 7: Make plots
%%
%Load in csv of sample data 
load('BMS01_sample_data')

%%
figure(2)
subplot(1,2,1)
semilogx(C_10Be_m,(z/100),'-b');
hold on
plot(C_10Be_m+C_10Be_std,(z/100),'--b');
plot(C_10Be_m-C_10Be_std,(z/100),'--b');
scatter(BMS01SampleData.C_10Be_ag,BMS01SampleData.depth_m,30,'k','filled');
eb_R(1) = errorbar(BMS01SampleData.C_10Be_ag,BMS01SampleData.depth_m,BMS01SampleData.C_10Be_unc_ag,...
    'horizontal', 'LineStyle', 'none');
eb_R(2) = errorbar(BMS01SampleData.C_10Be_ag,BMS01SampleData.depth_m,BMS01SampleData.depth_amalg, 'vertical', 'LineStyle', 'none');
set(eb_R, 'color', 'b', 'LineWidth', 1)
set(gca, 'YDir','reverse')
ylabel('Depth, m')
xlabel('^1^0Be, a/g')
ylim([0, 18.5])
subplot(1,2,2)
plot(Ratio_m,(z/100),'-k');
hold on
plot(Ratio_m+Ratio_std,(z/100),'--k');
plot(Ratio_m-Ratio_std,(z/100),'--k');
scatter(BMS01SampleData.Al10Be, BMS01SampleData.depth_m,30,'k','filled');
eb_R(1) = errorbar(BMS01SampleData.Al10Be,BMS01SampleData.depth_m,BMS01SampleData.Al10Be_unc,...
    'horizontal', 'LineStyle', 'none');
eb_R(2) = errorbar(BMS01SampleData.Al10Be,BMS01SampleData.depth_m,BMS01SampleData.depth_amalg, 'vertical', 'LineStyle', 'none');
set(eb_R, 'color', 'k', 'LineWidth', 1)
set(gca, 'YDir','reverse')
xlabel('^2^6Al / ^1^0Be')
ylim([0, 18.5])
