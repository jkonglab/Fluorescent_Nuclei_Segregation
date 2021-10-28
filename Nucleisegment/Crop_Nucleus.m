clc;
close all;

%I=imread('/Users/xiaoyuanguo/Desktop/Data/test_nuclei8.tif');
I=imread( '/Users/xiaoyuanguo/Desktop/7_OS15-915/All/OS15-915 63x_BT_2.tif');

% I=imread('/Users/xiaoyuanguo/Documents/Matlab_project/MATLAB/nucleiData/S16_28151 63x_PP_2.tif');
I2 = imcrop(I);
imshow(I2)
imwrite(I2,'/Users/xiaoyuanguo/Desktop/Data/test_nuclei52.tif');