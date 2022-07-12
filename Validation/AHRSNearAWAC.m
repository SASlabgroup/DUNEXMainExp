%% AHRS near AWAC Test
clear all
close all
clc

% Load the data from each mission
data_file = ["nearAWACData_mission18_microSWIFT_4.mat"];
data_1 = load(data_file);

%% Set up AHRS filter
fs = 12;
dt = 1/fs;
% **** NOTE that I changed it to NED from ENU since our sensors
% are measuring in that frame of ref
fuse = ahrsfilter('SampleRate', fs,'ReferenceFrame','NED','DecimationFactor',1);
fuse.AccelerometerNoise = 1e-4; % 1e-5 (m/s²)² 
fuse.GyroscopeNoise = deg2rad(5e0); % deg2rad(5e-2) (rad/s)²
fuse.MagnetometerNoise = 1; % 0.5 (µT)²
%fuse.GyroscopeDriftNoise = 1e-5; %(rad/s)²
%fuse.LinearAccelerationNoise = 1e-8; % (m/s²)² 
accel = transpose([data_1.a_x; data_1.a_y; data_1.a_z]);
gyro = transpose([data_1.gyro_x; data_1.gyro_y; data_1.gyro_z]);
mag = transpose([data_1.mag_x; data_1.mag_y; data_1.mag_z]);
[ q , rotv ] = fuse(accel, deg2rad(gyro), mag );

%% Compute Euler Angles from the Quaternions
ENU.angles = eulerd(q,'XYZ','point');
time = linspace(0,length(data_1.time), length(data_1.time));
ENU.time = time;

figure(1), clf, 
plot(ENU.time,ENU.angles,'.'), ylabel('Euler angles [deg]'), datetick

%% rotate accelerations
R = rotmat(q,'point');
for ri=1:length(time)
    ENU.acc(ri,:) = squeeze(R(:,:,ri)) * accel(ri,:)';
end

figure(2), clf
subplot(1,2,1), plot(ENU.acc),  ylabel('ENU acceleration [m/s^2]'),
subplot(1,2,2), pwelch(detrend(ENU.acc(round(60*fs):end,:)),[],[],[],fs), set(gca,'XScale','log')

%% Compare the rotated and original accelerations
figure(3)
% a_x
subplot(3,1,1)
plot(accel(:,1))
hold on
plot(ENU.acc(:,1))
ylabel('Acceleration [m/s^2]')
xlabel('Time Index')
title('East-West')
legend('original', 'corrected')

% a_y
subplot(3,1,2)
plot(accel(:,2))
hold on
plot(ENU.acc(:,2))
ylabel('Acceleration [m/s^2]')
xlabel('Time Index')
title('North-South')
legend('original', 'corrected')

% a_z
subplot(3,1,3)
plot(accel(:,3))
hold on
plot(ENU.acc(:,3))
ylabel('Acceleration [m/s^2]')
xlabel('Time Index')
title('Up-Down')
legend('original', 'corrected')

%% Different in verrtical Acceleration
accel_diff = ENU.acc(:,3) - accel(:,3);
figure(4)
plot(accel_diff)
ylabel('Acceleration Difference')
xlabel('Time Index')

%% Save the corrected data
ax_corrected = ENU.acc(:,1);
ay_corrected = ENU.acc(:,2);
az_corrected = ENU.acc(:,3);
ax_original = accel(:,1);
ay_original = accel(:,2);
az_original = accel(:,3);
save('ahrs_mission18_microSWIFT_4.mat', "ax_corrected", "ay_corrected", "az_corrected", "ax_original", "ay_original", "az_original")


