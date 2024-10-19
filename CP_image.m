

clear all;
close all;
clc;

image1 = imread('Versova Beach (1).png'); 
image2 = imread('Versova Beach (2).png'); 

image1_gray = rgb2gray(image1);
image2_gray = rgb2gray(image2);

edges1 = edge(image1_gray, 'Canny');
edges2 = edge(image2_gray, 'Canny');

difference = abs(double(edges2) - double(edges1));

threshold = 0.5; 
binary_mask = difference > threshold;

color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);

output = image1;
output(color_mask) = image2(color_mask);

figure;

imshow(output);
title('Output Image');
shoreline1 = sum(edges1(:));
shoreline2 = sum(edges2(:));

erosion_length_pixels = shoreline2 - shoreline1;

disp(['Erosion Length (Pixels): ' num2str(erosion_length_pixels)]);

conversion_factor = 0.1;
erosion_length_feet = erosion_length_pixels * conversion_factor;

disp(['Erosion Length (Feet): ' num2str(erosion_length_feet)]);

shoreline_change_percentage = (erosion_length_pixels / shoreline1) * 100;

disp(['Percentage of Shoreline Change: ' num2str(shoreline_change_percentage) '%']);

years = 1; 
average_retreat_rate_feet_per_year = erosion_length_feet / years;

disp(['Average Retreat Rate (Feet/Year): ' num2str(average_retreat_rate_feet_per_year)]);

combined = imfuse(image1, image2, 'blend', 'Scaling', 'joint');

difference = abs(double(edges2) - double(edges1));

threshold = graythresh(difference);
erosionMask = imbinarize(difference, threshold);

combinedErosion = combined;
combinedErosion(repmat(~erosionMask, 1, 1, 3)) = 0;

figure;
subplot(2, 2, 1);
imshow(image1);
title('Image 1');

subplot(2, 2, 2);
imshow(image2);
title('Image 2');

subplot(2, 2, 3);
imshow(difference);
title('Difference');

subplot(2, 2, 4);
imshow(combinedErosion);
title('Combined with Erosion');




