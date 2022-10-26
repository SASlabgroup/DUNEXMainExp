function [accel_x_earth, accel_y_earth, accel_z_earth] = AHRSAccelCorrection(accel_x_body, accel_y_body, accel_z_body, gyro_x_body, gyro_y_body, gyro_z_body, mag_x_body, mag_y_body, mag_z_body)
% AHRSAccelCorrection: Written by: EJ Rainville, Summer 2022
% Description: This function takes in measurements from an IMU onboard the microSWIFT and corrects the 
% accelerations into the earth frame of reference from the frame of reference of the buoy.

% Set up AHRS filter
fs = 12;
fuse = ahrsfilter('SampleRate', fs,'ReferenceFrame','NED','DecimationFactor',1);
fuse.AccelerometerNoise = 0.00004; % Noise variance from accelerometer, units are (m/s²)² 
fuse.GyroscopeNoise = deg2rad(0.045); % Noise variance from gyroscope, units are deg2rad(5e-2) (rad/s)²
fuse.MagnetometerNoise = 2; % Noise from Magnetometer, units are 0.5 (µT)²

% Organize data structures of accelerations, rotation rates and headings
accel = transpose([accel_x_body; accel_y_body; accel_z_body]);
gyro = transpose([gyro_x_body; gyro_y_body; gyro_z_body]);
mag = transpose([mag_x_body; mag_y_body; mag_z_body]);

% Fuse the IMU measurements and estimate the heading quaternion, q
[ q , ~ ] = fuse(accel, deg2rad(gyro), mag);

% Rotate Accelerations with euler angles 
R = rotmat(q,'point');
accel_earth = zeros(size(accel));
for ri=1:length(accel_x_body)
    accel_earth(ri,:) = squeeze(R(:,:,ri)) * accel(ri,:)';
end

% Define Rotated accelerations
accel_x_earth = accel_earth(:,1)';
accel_y_earth = accel_earth(:,2)';
accel_z_earth = accel_earth(:,3)';

end