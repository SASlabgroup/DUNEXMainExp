% load data
clear all; clc; 
data_file = ["./data/nearAWACData_mission19_microSWIFT_23.mat"];
data = load(data_file);

% organize variables
accel_x_body = data.a_x;
accel_y_body = data.a_y; 
accel_z_body = data.a_z; 
gyro_x_body = data.gyro_x; 
gyro_y_body = data.gyro_y; 
gyro_z_body = data.gyro_z;
mag_x_body = data.mag_x; 
mag_y_body = data.mag_y; 
mag_z_body = data.mag_z; 

% test function
[accel_x_earth, accel_y_earth, accel_z_earth] = AHRSAccelCorrection(accel_x_body, accel_y_body, accel_z_body, gyro_x_body, gyro_y_body, gyro_z_body, mag_x_body, mag_y_body, mag_z_body);

% Plot accels
figure
plot(accel_z_body, 'r')
hold on
plot(accel_z_earth, 'b')