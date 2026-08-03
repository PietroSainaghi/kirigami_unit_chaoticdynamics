function [ERR, ERR_train, ERR_test, recognition_train, recognition_test] = test_train_PhysicallyUnique(PARAMS, nGeo, nTrain, nTest, nExperiment, metric)




%% compute compound metric vectors

for geoIDX = 1:nGeo
    for mountIDX = 1:nTrain
        for metricIDX = 1:size(metric,3)
            metricVAL_train(geoIDX,mountIDX,metricIDX) = mean( metric(geoIDX, ((mountIDX-1)*nExperiment+1):(nExperiment*mountIDX),metricIDX) );
            metricSTD_train(geoIDX,mountIDX,metricIDX) = std( metric(geoIDX, ((mountIDX-1)*nExperiment+1):(nExperiment*mountIDX),metricIDX) );
        end
    end
end



for geoIDX = 1:nGeo
    for mountIDX = 1:nTest
        for metricIDX = 1:size(metric,3)
            metricVAL_test(geoIDX,mountIDX,metricIDX) = mean( metric(geoIDX, ((mountIDX-1)*nExperiment+1+nExperiment*nTrain):(nExperiment*mountIDX+nExperiment*nTrain) ,metricIDX) );
            metricSTD_test(geoIDX,mountIDX,metricIDX) = std( metric(geoIDX, ((mountIDX-1)*nExperiment+1+nExperiment*nTrain):(nExperiment*mountIDX+nExperiment*nTrain) ,metricIDX) );
        end
    end
end

%% set threshold

for geoIDX = 1:nGeo
    for metricIDX = 1:size(metric,3)
        metricSTD_geo = sqrt( sum( ( metricSTD_train(geoIDX,:,metricIDX) ).^2 ) );
        threshold(geoIDX,metricIDX) = PARAMS( (geoIDX-1)*size(metric,3) + metricIDX ) .* metricSTD_geo;
    end
end


%% compute metric error

for geoIDX1 = 1:nGeo
    for mountIDX1 = 1:nTrain
        for geoIDX2 = 1:nGeo
            for mountIDX2 = 1:nTrain
                for metricIDX = 1:size(metric,3)
                    SE_metric_train( (geoIDX1-1)*nTrain+mountIDX1, (geoIDX2-1)*nTrain+mountIDX2, metricIDX ) = abs( metricVAL_train(geoIDX1,mountIDX1,metricIDX) - metricVAL_train(geoIDX2,mountIDX2,metricIDX) );
                end
            end
        end
    end
end


for geoIDX1 = 1:nGeo
    for mountIDX1 = 1:nTrain
        for geoIDX2 = 1:nGeo
            for mountIDX2 = 1:nTest
                for metricIDX = 1:size(metric,3)
                    SE_metric_test( (geoIDX1-1)*nTrain+mountIDX1, (geoIDX2-1)*nTest+mountIDX2, metricIDX ) = abs( metricVAL_train(geoIDX1,mountIDX1,metricIDX) - metricVAL_test(geoIDX2,mountIDX2,metricIDX) );
                end
            end
        end
    end
end

%% compute pass/fail

k_sigmoid = 200;

for geomountIDX1 = 1:size(SE_metric_train,1)
    for geomountIDX2 = 1:size(SE_metric_train,2)
        for metricIDX = 1:size(SE_metric_train,3)
            % passFail_train(geomountIDX1,geomountIDX2,metricIDX) =  SE_metric_train(geomountIDX1,geomountIDX2,metricIDX) < threshold(ceil(geomountIDX1/nTrain),metricIDX);
            diff_val_train = SE_metric_train(geomountIDX1,geomountIDX2,metricIDX) - threshold(ceil(geomountIDX1/nTrain),metricIDX);
            passFail_train(geomountIDX1,geomountIDX2,metricIDX) = 1 / (1 + exp(k_sigmoid * diff_val_train));
        end
    end
end

for geomountIDX1 = 1:size(SE_metric_test,1)
    for geomountIDX2 = 1:size(SE_metric_test,2)
        for metricIDX = 1:size(SE_metric_test,3)
            % passFail_test(geomountIDX1,geomountIDX2,metricIDX) = SE_metric_test(geomountIDX1,geomountIDX2,metricIDX) < threshold(ceil(geomountIDX1/nTrain),metricIDX);
            diff_val_test = SE_metric_test(geomountIDX1,geomountIDX2,metricIDX) - threshold(ceil(geomountIDX1/nTrain),metricIDX);
            passFail_test(geomountIDX1,geomountIDX2,metricIDX) = 1 / (1 + exp(k_sigmoid * diff_val_test));
        end
    end
end


%% compute recognition matrix

for geoIDX1 = 1:nGeo
    for geoIDX2 = 1:nGeo
        recognition_train(geoIDX1,geoIDX2) = sum(sum(sum( passFail_train( ((geoIDX1-1)*nTrain+1):(geoIDX1*nTrain) , ((geoIDX2-1)*nTrain+1):(geoIDX2*nTrain) , : ) ))) / (nTrain*nTrain*size(metric,3));
        recognition_test(geoIDX1,geoIDX2) = sum(sum(sum( passFail_test( ((geoIDX1-1)*nTrain+1):(geoIDX1*nTrain) , ((geoIDX2-1)*nTest+1):(geoIDX2*nTest) , : ) ))) / (nTrain*nTest*size(metric,3));
    end
end

%% compute error

for geoIDX1 = 1:nGeo
    for geoIDX2 = 1:nGeo
        if geoIDX1 == geoIDX2
            diagScore_train(geoIDX1) = recognition_train(geoIDX1,geoIDX2);
            diagScore_test(geoIDX1) = recognition_test(geoIDX1,geoIDX2);
            offScore_train(geoIDX1,geoIDX2) = 0;
            offScore_test(geoIDX1,geoIDX2) = 0;
        else
            offScore_train(geoIDX1,geoIDX2) = recognition_train(geoIDX1,geoIDX2);
            offScore_test(geoIDX1,geoIDX2) = recognition_test(geoIDX1,geoIDX2);
        end
    end
end


%% iteration 1: reciprocal

% better than iteration 2 in test

%RegOff = 0.6;
%ERR_train = (1-RegOff)*mean( 1 ./ (diagScore_train.^2) ) + RegOff*mean(mean((offScore_train).^2));
%ERR_test = (1-RegOff)*mean( 1 ./ (diagScore_test.^2) ) + RegOff*mean(mean((offScore_test).^2));


%% iteration 2: difference

% great in testing, does overfit

% RegOff = 0.5;
% 
% ERR_train = (1-RegOff)*mean( (1 - diagScore_train.^2) ) + RegOff*mean(mean((offScore_train).^2));
% ERR_test = (1-RegOff)*mean( (1 - diagScore_test.^2) ) + RegOff*mean(mean((offScore_test).^2));
% 

%% iteration 3: margin loss

% 

margin = 0.3; 

genuine_train = diagScore_train;
imposter_train = offScore_train(~logical(eye(nGeo)));

loss_genuine_train = mean((1 - genuine_train).^2);
loss_imposter_train = mean(max(0, margin - (genuine_train - imposter_train)).^2, 'all');

genuine_test = diagScore_test;
imposter_test = offScore_test(~logical(eye(nGeo)));

loss_genuine_test = mean((1 - genuine_test).^2);
loss_imposter_test = mean(max(0, margin - (genuine_test - imposter_test)).^2, 'all');


RegOff = 0.4;

ERR_train = (1-RegOff)*loss_genuine_train + RegOff*loss_imposter_train;
ERR_test = (1-RegOff)*loss_genuine_test + RegOff*loss_imposter_test;





%% assignment

ERR = ERR_train;










