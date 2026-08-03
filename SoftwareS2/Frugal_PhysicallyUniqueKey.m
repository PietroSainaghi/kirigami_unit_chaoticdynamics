clear
clc
close all

%% load dataset

nGeo = 5;
nTrain = 3;
nTest = 2;
nExperiment = 10;

addpath('Data')

allFiles = dir('**');

Ngood = 0;
for fileIDX = 1:size(allFiles,1)
    if contains(allFiles(fileIDX,1).name, '.mat')
        if contains(allFiles(fileIDX,1).name, 'Key')
            Ngood = Ngood+1;
            matlabFiles(Ngood,1) = allFiles(fileIDX,1);
        end
    end
end

for fileIDX = 1:size(matlabFiles,1)
    geoIDX = ceil(fileIDX/((nTrain+nTest)*nExperiment));
    runIDX = mod(fileIDX,((nTrain+nTest)*nExperiment));
    if runIDX == 0
        runIDX = (nTrain+nTest)*nExperiment;
    end
    rawDataSet{fileIDX} = load([matlabFiles(fileIDX,1).folder,'\',matlabFiles(fileIDX,1).name]);
    t_u{geoIDX,runIDX} = rawDataSet{fileIDX}.moku.data(:,1);
    V_F_u{geoIDX,runIDX} = rawDataSet{fileIDX}.moku.data(:,2);
    V_v_u{geoIDX,runIDX} = rawDataSet{fileIDX}.moku.data(:,3);
end



%% convert from voltage into reading

% vibrometer gain
vG = 20; % m/(s*V)

% moku gain
mG = 10; % -20 dB correction

% load sensor gain
FG = 21.03e-3; % V/N


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        F_s{geoIDX,runIDX} = V_F_u{geoIDX,runIDX} / FG;
        v_s{geoIDX,runIDX} = V_v_u{geoIDX,runIDX} * vG * mG;
    end
end


%% get steady state

ss_start = 6;
ss_end = 10;

for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        [~,startIDX] = max((abs(t_u{geoIDX,runIDX}-ss_start)<1e-6));
        [~,endIDX] = max((abs(t_u{geoIDX,runIDX}-ss_end)<1e-6));
        v_ss{geoIDX,runIDX} = v_s{geoIDX,runIDX}(startIDX:endIDX);
        t_ss{geoIDX,runIDX} = t_u{geoIDX,runIDX}(startIDX:endIDX);
    end
end


%% integrate into position


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        z_ss{geoIDX,runIDX} = 1000*cumtrapz(t_ss{geoIDX,runIDX}, highpass(v_ss{geoIDX,runIDX},70,50000) );
    end
end


%% convert phase portrait to polar coordinates


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        v_smooth = smooth(v_ss{geoIDX,runIDX},15,'sgolay',9);
        z_smooth = smooth(z_ss{geoIDX,runIDX},15,'sgolay',9);
        theta_ss{geoIDX,runIDX} = atan2(v_smooth, z_smooth);
        r_ss{geoIDX,runIDX} = sqrt( (v_smooth).^2 + (z_smooth).^2 );
    end
end

%% fourier transform

Fs = 50000;


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        L = length(v_s{geoIDX,runIDX});
        fft_ss{geoIDX,runIDX} = abs(fftshift(fft(v_s{geoIDX,runIDX})))/max(abs(fftshift(fft(v_s{geoIDX,runIDX}))));
        f_ss{geoIDX,runIDX} = Fs/L*(-L/2:L/2-1);
    end
end


%% Metrics 1, 2, 3: max min mean velocity


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        metric(geoIDX,runIDX,1) = mean(findpeaks(v_ss{geoIDX,runIDX}));
        metric(geoIDX,runIDX,2) = mean(findpeaks(-v_ss{geoIDX,runIDX}));
        metric(geoIDX,runIDX,3) = mean(v_ss{geoIDX,runIDX});
    end
end

%% Metrics 4, 5, 6: max min mean displacement


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        metric(geoIDX,runIDX,4) = mean(findpeaks(z_ss{geoIDX,runIDX}));
        metric(geoIDX,runIDX,5) = mean(findpeaks(-z_ss{geoIDX,runIDX}));
        metric(geoIDX,runIDX,6) = mean(z_ss{geoIDX,runIDX});
    end
end


%% Metrics 7 - 15: harmonics in fourier domain

fourierViaCoords = 1850 .* [0, 1/16, 1/8, 1/4, 1/2, 1, 2, 3, 4];

for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        for viaIDX = 1:length(fourierViaCoords)
            df_fft = min(f_ss{geoIDX,runIDX}(2:end)-f_ss{geoIDX,runIDX}(1:end-1));
            viaIndices =  abs( f_ss{geoIDX,runIDX} - fourierViaCoords(viaIDX) ) <= (df_fft);
            if sum(viaIndices) == 1
                metric(geoIDX,runIDX,6+viaIDX) = fft_ss{geoIDX,runIDX}(viaIndices);
            elseif sum(viaIndices) == 0
                multiVia = fft_ss{geoIDX,6+runIDX}(abs( f_ss{geoIDX,runIDX} - fourierViaCoords(viaIDX) ) <= (df_fft*2));
                metric(geoIDX,runIDX,viaIDX) = multiVia(1);
            elseif sum(viaIndices) > 1
                multiVia = fft_ss{geoIDX,runIDX}(viaIndices);
                metric(geoIDX,runIDX,6+viaIDX) = multiVia(1);
            end
        end
    end
end


%% metric 16 17: chaos metrics for velocity


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        v_temp = v_ss{geoIDX,runIDX}((end-30000):end);
        [~, delay, dim] = phaseSpaceReconstruction(v_temp);
        if isempty(delay)
            delay = 1
        end
        if isempty(dim)
            dim = 4
        end
        metric(geoIDX, runIDX, 16) = lyapunovExponent(v_temp, Fs, delay, dim);
        metric(geoIDX, runIDX, 17) = correlationDimension(v_temp, delay, dim);
    end
end


%% metric 18, 19: chaos metrics for position


for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        z_temp = z_ss{geoIDX,runIDX}((end-30000):end);
        [~, delay, dim] = phaseSpaceReconstruction(z_temp);
        if isempty(delay)
            delay = 1
        end
        if isempty(dim)
            dim = 4
        end
        metric(geoIDX, runIDX, 18) = lyapunovExponent(z_temp, Fs, delay, dim);
        metric(geoIDX, runIDX, 19) = correlationDimension(z_temp, delay, dim);
    end
end

%% metric 20-22: dispersion for velocity

for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        metric(geoIDX, runIDX, 20) = rms(v_ss{geoIDX,runIDX});
        metric(geoIDX, runIDX, 21) = skewness(v_ss{geoIDX,runIDX});
        metric(geoIDX, runIDX, 22) = kurtosis(v_ss{geoIDX,runIDX});
    end
end

%% metric 23-25: dispersion for velocity

for geoIDX = 1:size(t_u,1)
    for runIDX = 1:size(t_u,2)
        metric(geoIDX, runIDX, 23) = rms(z_ss{geoIDX,runIDX});
        metric(geoIDX, runIDX, 24) = skewness(z_ss{geoIDX,runIDX});
        metric(geoIDX, runIDX, 25) = kurtosis(z_ss{geoIDX,runIDX});
    end
end



%% train model

options = optimoptions('ga',...
    'PlotFcn', @gaplotbestf,...
    'PopulationSize',500, ...
    'MaxGenerations', 10000,...
    'ConstraintTolerance',1e-20,...
    'FunctionTolerance',1e-20,...
    'MaxStallGenerations', 1000,...
    'UseParallel', true);

lb = 0.1.*ones(1,size(metric,3)*nGeo);
ub = 10.*ones(1,size(metric,3)*nGeo);


[x_params,fval,exitflag,output,population,scores] = ga( ...
        @(PARAMS) test_train_PhysicallyUnique(PARAMS, nGeo, nTrain, nTest, nExperiment, metric),...
        length(ub),...
        [],...
        [],...
        [],...
        [],...
        lb,...
        ub,...
        [],...
        [],...
        options);


[ERR, ERR_train, ERR_test, recognition_train, recognition_test] = test_train_PhysicallyUnique(x_params, nGeo, nTrain, nTest, nExperiment, metric);


recognition_train


recognition_test




