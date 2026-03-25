close all;
clear;

data = readtable('ibtracs.NA.list.v04r01 (1).csv');

data(strcmp(data.NAME, 'UNNAMED'), :) = [];

long = data.LON;
lat = data.LAT;

figure();
plot(long(1:38), lat(1:38))

% storms = data(strcmp(data.NAME, "storms") & data.SEASON == 2012, :);
storms = data(data.SEASON==2000, :);

pred_lon = zeros(size(storms, 1), 1);
pred_lat = zeros(size(storms, 1), 1);

pred_lon(1) = storms.LON(1);
pred_lat(1) = storms.LAT(1);

pred_lon(2) = storms.LON(2);
pred_lat(2) = storms.LAT(2);


pred_lon_com(1:2) = pred_lon(1:2);
pred_lat_com(1:2) = pred_lon(1:2);
error_lon(1:2) = 0;
error_lat(1:2) = 0;

dt = 3; % Timestep in hours (since speed is mph)

idx = 3;

storm_ids = strcat(storms.NAME, '_', string(storms.SEASON));

storm_size = length(storm_ids);

for storm=1:storm_size
    for i=2:size(storms,1)-1
        pred_lon(i+1) = storms.LON(i) + dt * (storms.LON(i)-storms.LON(i-1));
        pred_lat(i+1) = storms.LAT(i) + dt * (storms.LAT(i)-storms.LAT(i-1));
        
        pred_lon_com(idx) = pred_lon(i+1);
        pred_lat_com(idx) = pred_lat(i+1);
        error_lon(idx) = (storms.LON(i+1) - pred_lon(i+1));
        error_lat(idx) = (storms.LAT(i+1) - pred_lat(i+1));
        
        idx = idx + 1;
    end
end

X = [pred_lon_com.', pred_lat_com.', ones(idx-1, 1)];

lon_coef = X \ ((error_lon.').^2);
lat_coef = X \ ((error_lat.').^2);

figure()
plot(storms.LON, storms.LAT, 'b');
hold on;
plot(pred_lon, pred_lat, 'r')
legend("Real Path", "First Order Prediction")
xlabel("Longitude")
ylabel("Latitude")
title("Predicting Hurricane Sandy's Path")
grid on

% hold on;
% plot(pred_lon + (error_lon.').^2, pred_lat + (error_lat.').^2, '--')
% hold on
% plot(pred_lon - (error_lon.').^2, pred_lat - (error_lat.').^2, '--')

% GIVEN THE FIRST 5 ENTRIES, PREDICT SANDY

sandy = data(strcmp(data.NAME, "SANDY") & data.SEASON == 2012, :);
observed_lon = sandy.LON(1:6);
observed_lat = sandy.LAT(1:6);

steps = 56; % 1 Week

sandy_pred_lat = zeros(steps);
sandy_pred_lon = zeros(steps);
sandy_err_lat = zeros(steps);
sandy_err_lon = zeros(steps);

for i=1:steps
    lon_cur = observed_lon(end);
    lat_cur = observed_lat(end);
    lon_last  = observed_lon(end-1);
    lat_last  = observed_lat(end-1);

    sandy_pred_lon(i) = lon_cur + dt * (lon_cur - lon_last);
    sandy_pred_lat(i) = lat_cur + dt * (lat_cur - lat_last);

    sandy_err_lat(i) = lat_coef(1) * lon_cur + lat_coef(2) * lat_cur + lat_coef(3);
    sandy_err_lon(i) = lon_coef(1) * lon_cur + lon_coef(2) * lat_cur + lon_coef(3);

    observed_lon = [observed_lon; sandy_pred_lon(i)];
    observed_lat = [observed_lat; sandy_pred_lat(i)];
   
end

figure()
plot(sandy_err_lon, sandy_err_lat, '.');