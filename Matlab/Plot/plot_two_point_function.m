clear; clc;
close all;

filename = '../../cmake-build-debug/two_point_functions';
file_id = fopen(filename,'rb');

row_data = fread(file_id,'double');

fclose(file_id);