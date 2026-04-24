close all;
clear;

data = readtable('ibtracs.NA.list.v04r01.csv');

data(strcmp(data.NAME, 'UNNAMED'), :) = [];

long = data.LON;
lat = data.LAT;

figure();
plot(long(1:38), lat(1:38))

% storms = data(strcmp(data.NAME, "storms") & data.SEASON == 2012, :);
storms = data(data.SEASON>=1990, :);

pred_lon = zeros(size(storms, 1), 1);
pred_lat = zeros(size(storms, 1), 1);

pred_lon(1) = storms.LON(1);
pred_lat(1) = storms.LAT(1);

pred_lon(2) = storms.LON(2);
pred_lat(2) = storms.LAT(2);

% Feature at each training sample is X_t = (lon_t, lat_t). We model
% log(squared error) as a linear function of X_t so that the implied
% variance exp(A*lon + B*lat + C) is always positive.
feat_lon = [];
feat_lat = [];
error_lon = [];
error_lat = [];

dt = 1;    % Timestep in hours (since speed is mph)
H  = 8;    % Forecast horizon (in steps) used to generate training targets.
           % IBTrACS positions are rounded to 0.1 deg, so one-step errors
           % frequently hit the quantization floor. H-step errors are big
           % enough to be resolved cleanly.

storm_ids = strcat(storms.NAME, '_', string(storms.SEASON));
unique_storm_ids = unique(storm_ids, 'stable');

for s = 1:length(unique_storm_ids)
    storm_mask = strcmp(storm_ids, unique_storm_ids{s});
    storm_lons = storms.LON(storm_mask);
    storm_lats = storms.LAT(storm_mask);
    storm_indices = find(storm_mask);

    if length(storm_lons) < 3
        continue;
    end

    for k = 2:length(storm_lons)-1
        i = storm_indices(k);
        % One-step prediction (kept for plotting, not used in training)
        pred_lon(i+1) = storm_lons(k) + dt * (storm_lons(k)-storm_lons(k-1));
        pred_lat(i+1) = storm_lats(k) + dt * (storm_lats(k)-storm_lats(k-1));

        % H-step training sample (only when there is enough lookahead).
        % Constant-velocity free-run: after H steps the prediction is
        %   X_k + H * (X_k - X_{k-1})
        if k + H <= length(storm_lons)
            v_lon = storm_lons(k) - storm_lons(k-1);
            v_lat = storm_lats(k) - storm_lats(k-1);
            pred_lon_H = storm_lons(k) + H * v_lon;
            pred_lat_H = storm_lats(k) + H * v_lat;

            feat_lon(end+1)  = storm_lons(k);
            feat_lat(end+1)  = storm_lats(k);
            error_lon(end+1) = storm_lons(k+H) - pred_lon_H;
            error_lat(end+1) = storm_lats(k+H) - pred_lat_H;
        end
    end
end

% Fit log(H-step variance) as a linear function of X_t:
% log(var_H) = A*lon + B*lat + C   ->   var_H = exp(A*lon + B*lat + C)
eps_safe = 1e-6;
X = [feat_lon.', feat_lat.', ones(length(feat_lon), 1)];

lon_coef = X \ log((error_lon.').^2 + eps_safe);
lat_coef = X \ log((error_lat.').^2 + eps_safe);

figure()
sandy = data(strcmp(data.NAME, "SANDY") & data.SEASON == 2012, :);
plot(sandy.LON, sandy.LAT, 'b');
hold on;

legend("Real Path", "First Order Prediction")
xlabel("Longitude")
ylabel("Latitude")
title("Real World Hurricane Sandy's Path")
grid on

% GIVEN THE FIRST 5 ENTRIES, PREDICT SANDY

sandy = data(strcmp(data.NAME, "SANDY") & data.SEASON == 2012, :);
observed_lon = sandy.LON(1:6);
observed_lat = sandy.LAT(1:6);

steps = 56; % 1 Week

sandy_pred_lat = zeros(steps, 1);
sandy_pred_lon = zeros(steps, 1);
sandy_err_lat = zeros(steps, 1);
sandy_err_lon = zeros(steps, 1);

for i=1:steps
    lon_cur = observed_lon(end);
    lat_cur = observed_lat(end);
    lon_last  = observed_lon(end-1);
    lat_last  = observed_lat(end-1);

    sandy_pred_lon(i) = lon_cur + dt * (lon_cur - lon_last);
    sandy_pred_lat(i) = lat_cur + dt * (lat_cur - lat_last);

    % Model gives H-step-ahead variance at the current position.
    % Under the independent-step assumption, per-step variance is
    % H-step variance divided by H.
    var_H_lon = exp(lon_coef(1) * lon_cur + lon_coef(2) * lat_cur + lon_coef(3));
    var_H_lat = exp(lat_coef(1) * lon_cur + lat_coef(2) * lat_cur + lat_coef(3));
    sandy_err_lon(i) = var_H_lon / H;
    sandy_err_lat(i) = var_H_lat / H;

    observed_lon = [observed_lon; sandy_pred_lon(i)];
    observed_lat = [observed_lat; sandy_pred_lat(i)];
   
end

% Propagate per-step variance forward (independent-step assumption).
% Cumulative variance at step n is the sum of per-step variances 1..n,
% so the uncertainty band grows with forecast horizon.
sandy_cum_var_lon = cumsum(sandy_err_lon);
sandy_cum_var_lat = cumsum(sandy_err_lat);

figure()
plot(sqrt(sandy_cum_var_lon), sqrt(sandy_cum_var_lat), '.');

figure()
plot(sandy_pred_lon + sqrt(sandy_cum_var_lon), sandy_pred_lat + sqrt(sandy_cum_var_lat), 'b');
hold on
plot(sandy_pred_lon - sqrt(sandy_cum_var_lon), sandy_pred_lat - sqrt(sandy_cum_var_lat), 'b');
hold on
plot(sandy_pred_lon, sandy_pred_lat)
title("Predicted Path of Huricane Sandy")

figure()
plot(sandy_pred_lon, sandy_pred_lat, 'k', 'LineWidth', 1.5)
hold on

% Draw an ellipse at every Nth step (otherwise the plot gets cluttered)
theta = linspace(0, 2*pi, 60);
ellipse_stride = 8;
for i = 1:ellipse_stride:steps
    a = sqrt(sandy_cum_var_lon(i));  % semi-axis in longitude
    b = sqrt(sandy_cum_var_lat(i));  % semi-axis in latitude
    ex = sandy_pred_lon(i) + a * cos(theta);
    ey = sandy_pred_lat(i) + b * sin(theta);
    plot(ex, ey, 'Color', [0.2 0.4 0.8 0.3]);  % translucent blue
end

% Outer envelope: connect the lat/lon extremes along the track
plot(sandy_pred_lon + sqrt(sandy_cum_var_lon), sandy_pred_lat, 'b');
plot(sandy_pred_lon - sqrt(sandy_cum_var_lon), sandy_pred_lat, 'b');
plot(sandy_pred_lon, sandy_pred_lat + sqrt(sandy_cum_var_lat), 'g');
plot(sandy_pred_lon, sandy_pred_lat - sqrt(sandy_cum_var_lat), 'g');

xlabel("Longitude")
ylabel("Latitude")
title("Sandy Predicted Path w/ Uncertainty")
grid on
