% SL1 = imread('Screenshot (623).png');
% SL2 = imread('Screenshot (624).png');
% SL12 = rgb2gray(SL1);
% SL22 = rgb2gray(SL2);
% RGB543 = cat(2, SL22, SL12);
% figure; imshow(RGB543);
% SL12 = double(SL12);
% SL22 = double(SL22);
% ndvi = SL22/SL12;
% figure; imshow(ndvi, []);
% colormap("jet");
% colorbar;
% [w, h] = size(SL12);
% acc1 = input("Please enter the threshold value = ");
% for i = 1:w
%     for j = 1:h
%         if ndvi(i,j)>=acc1
%             Smap(i,j,1) = 65535;
%             Smap(i,j,2) = 0;
%             Smap(i,j,3) = 0;
% 
%         else
%             Smap(i,j,1) = SL12(i,j);
%             Smap(i,j,2) = SL12(i,j);
%             Smap(i,j,3) = SL12(i,j);
%         end
%     end
% end
% Smap = uint16(Smap);
% figure; imshow(Smap);

% % Read the two images
% image2 = imread('Screenshot (623).png'); % Replace with the actual path to your first image
% image1 = imread('Screenshot (624).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Calculate the difference between the images
% erosion = abs(double(image2_gray) - double(image1_gray));
% 
% % Apply edge detection
% edges = edge(erosion, 'Canny');
% 
% % Apply thresholding to identify areas of erosion
% threshold = graythresh(erosion);
% erosionMask = imbinarize(erosion, threshold);
% 
% % Display the erosion results
% figure;
% subplot(1, 3, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(1, 3, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(1, 3, 3);
% imshow(erosionMask);
% title('Erosion Mask');

% -------------------------------------------------------------------------------------------------
% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Screenshot (623).png'); % Replace with the actual path to your first image
% image2 = imread('Screenshot (624).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Combine the two images
% combined = imfuse(image1, image2, 'blend', 'Scaling', 'joint');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Apply thresholding to identify areas of erosion
% threshold = graythresh(difference);
% erosionMask = imbinarize(difference, threshold);
% 
% % Apply the erosion mask to the combined image
% combinedErosion = combined;
% combinedErosion(repmat(~erosionMask, 1, 1, 3)) = 0;
% 
% % Display the results
% figure;
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
% imshow(combinedErosion);
% title('Combined with Erosion');
% -----------------------------------------------------------------------------------------

% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Screenshot (623).png'); % Replace with the actual path to your first image
% image2 = imread('Screenshot (624).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Threshold the difference image to obtain the binary mask
% threshold = 0.5; % Adjust the threshold as needed
% binary_mask = difference > threshold;
% 
% % Create a color mask based on the binary mask
% color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);
% 
% % Combine the color mask with the original image1
% output = image1;
% output(color_mask) = image2(color_mask);
% 
% % Display the output image
% figure;
% imshow(output);
% title('Output Image');
% ------------------------------------------------------------------------------------------------

% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Screenshot (627).png'); % Replace with the actual path to your first image
% image2 = imread('Screenshot (628).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Threshold the difference image to obtain the binary mask
% threshold = 0.5; % Adjust the threshold as needed
% binary_mask = difference > threshold;
% 
% % Create a color mask based on the binary mask
% color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);
% 
% % Combine the color mask with the original image1
% output = image1;
% output(color_mask) = image2(color_mask);
% 
% % Display the output image
% figure;
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
% imshow(output);
% title('Output Image');
% 
% % Combine the two images
% combined = imfuse(image1, image2, 'blend', 'Scaling', 'joint');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Apply thresholding to identify areas of erosion
% threshold = graythresh(difference);
% erosionMask = imbinarize(difference, threshold);
% 
% % Apply the erosion mask to the combined image
% combinedErosion = combined;
% combinedErosion(repmat(~erosionMask, 1, 1, 3)) = 0;
% 
% % Display the results
% figure;
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
% imshow(combinedErosion);
% title('Combined with Erosion');

% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------

% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Om Beach (1).png'); % Replace with the actual path to your first image
% image2 = imread('Om Beach (2).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Threshold the difference image to obtain the binary mask
% threshold = 0.5; % Adjust the threshold as needed
% binary_mask = difference > threshold;
% 
% % Create a color mask based on the binary mask
% color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);
% 
% % Combine the color mask with the original image1
% output = image1;
% output(color_mask) = image2(color_mask);
% 
% % Display the output image
% figure;
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
% imshow(output);
% title('Output Image');
% 
% % Calculate the length of the shoreline in both images
% shoreline1 = sum(edges1(:));
% shoreline2 = sum(edges2(:));
% 
% % Calculate the erosion length in pixels
% erosion_length_pixels = shoreline2 - shoreline1;
% 
% % Display the erosion length in pixels
% disp(['Erosion Length (Pixels): ' num2str(erosion_length_pixels)]);
% 
% % Convert erosion length to the desired scale (e.g., feet)
% % Assuming a conversion factor of 1 pixel = 0.1 feet
% conversion_factor = 0.1;
% erosion_length_feet = erosion_length_pixels * conversion_factor;
% 
% % Display the erosion length in feet
% disp(['Erosion Length (Feet): ' num2str(erosion_length_feet)]);
% 
% % Combine the two images
% combined = imfuse(image1, image2, 'blend', 'Scaling', 'joint');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Apply thresholding to identify areas of erosion
% threshold = graythresh(difference);
% erosionMask = imbinarize(difference, threshold);
% 
% % Apply the erosion mask to the combined image
% combinedErosion = combined;
% combinedErosion(repmat(~erosionMask, 1, 1, 3)) = 0;
% 
% % Display the results
% figure;
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
% imshow(combinedErosion);
% title('Combined with Erosion');

% ------------------------------------------------
% ------------------------------------------------
% ------------------------------------------------

% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Kappad Beach (1).png'); % Replace with the actual path to your first image
% image2 = imread('Kappad Beach (2).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Threshold the difference image to obtain the binary mask
% threshold = 0.5; % Adjust the threshold as needed
% binary_mask = difference > threshold;
% 
% % Create a color mask based on the binary mask
% color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);
% 
% % Combine the color mask with the original image1
% output = image1;
% output(color_mask) = image2(color_mask);
% 
% % Display the output image
% figure;
% imshow(output);
% title('Output Image');
% 
% % Calculate the length of the shoreline in both images
% shoreline1 = sum(edges1(:));
% shoreline2 = sum(edges2(:));
% 
% % Calculate the erosion length in pixels
% erosion_length_pixels = shoreline2 - shoreline1;
% 
% % Display the erosion length in pixels
% disp(['Erosion Length (Pixels): ' num2str(erosion_length_pixels)]);
% 
% % Convert erosion length to the desired scale (e.g., feet)
% % Assuming a conversion factor of 1 pixel = 0.1 feet
% conversion_factor = 0.1;
% erosion_length_feet = erosion_length_pixels * conversion_factor;
% 
% % Display the erosion length in feet
% disp(['Erosion Length (Feet): ' num2str(erosion_length_feet)]);
% 
% % Calculate the erosion area in pixels
% erosion_area_pixels = sum(binary_mask(:));
% 
% % Convert erosion area to the desired scale (e.g., square meters)
% % Assuming a conversion factor of 1 pixel = 0.25 square meters
% conversion_factor_area = 0.25;
% erosion_area_meters = erosion_area_pixels * conversion_factor_area;
% 
% % Display the erosion area in square meters
% disp(['Erosion Area (Square Meters): ' num2str(erosion_area_meters)]);

% 
% 
% 
%

% To compare more than 2 images

% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Screenshot (644).png'); % Replace with the actual path to your first image
% image2 = imread('Screenshot (645).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Threshold the difference image to obtain the binary mask
% threshold = 0.5; % Adjust the threshold as needed
% binary_mask = difference > threshold;
% 
% % Create a color mask based on the binary mask
% color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);
% 
% % Combine the color mask with the original image1
% output = image1;
% output(color_mask) = image2(color_mask);
% 
% % Display the output image
% figure;
% imshow(output);
% title('Output Image');
% 
% % Calculate the length of the shoreline in both images
% shoreline1 = sum(edges1(:));
% shoreline2 = sum(edges2(:));
% 
% % Calculate the erosion length in pixels
% erosion_length_pixels = shoreline2 - shoreline1;
% 
% % Display the erosion length in pixels
% disp(['Erosion Length (Pixels): ' num2str(erosion_length_pixels)]);
% 
% % Convert erosion length to the desired scale (e.g., feet)
% % Assuming a conversion factor of 1 pixel = 0.1 feet
% conversion_factor = 0.1;
% erosion_length_feet = erosion_length_pixels * conversion_factor;
% 
% % Display the erosion length in feet
% disp(['Erosion Length (Feet): ' num2str(erosion_length_feet)]);
% 
% % Calculate the erosion area in pixels
% erosion_area_pixels = sum(binary_mask(:));
% 
% % Convert erosion area to the desired scale (e.g., square meters)
% % Assuming a conversion factor of 1 pixel = 0.25 square meters
% conversion_factor_area = 0.25;
% erosion_area_meters = erosion_area_pixels * conversion_factor_area;
% 
% % Display the erosion area in square meters
% disp(['Erosion Area (Square Meters): ' num2str(erosion_area_meters)]);
% 
% % Calculate the average retreat rate of the shoreline in feet per year
% time_interval = 1; % Time interval between the two images in years
% retreat_rate_feet_per_year = erosion_length_feet / time_interval;
% 
% % Display the average retreat rate
% disp(['Average Retreat Rate (Feet/Year): ' num2str(retreat_rate_feet_per_year)]);


% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Aksa Beach (1).png'); % Replace with the actual path to your first image
% image2 = imread('Aksa Beach (2).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Threshold the difference image to obtain the binary mask
% threshold = 0.5; % Adjust the threshold as needed
% binary_mask = difference > threshold;
% 
% % Create a color mask based on the binary mask
% color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);
% 
% % Combine the color mask with the original image1
% output = image1;
% output(color_mask) = image2(color_mask);
% 
% % Display the output image
% figure;
% imshow(output);
% title('Output Image');
% 
% % Calculate the length of the shoreline in both images
% shoreline1 = sum(edges1(:));
% shoreline2 = sum(edges2(:));
% 
% % Calculate the erosion length in pixels
% erosion_length_pixels = shoreline2 - shoreline1;
% 
% % Display the erosion length in pixels
% disp(['Erosion Length (Pixels): ' num2str(erosion_length_pixels)]);
% 
% % Convert erosion length to the desired scale (e.g., meters)
% % Assuming a conversion factor of 1 pixel = 0.1 meters
% conversion_factor = 0.1;
% erosion_length_meters = erosion_length_pixels * conversion_factor;
% 
% % Display the erosion length in meters
% disp(['Erosion Length (Meters): ' num2str(erosion_length_meters)]);
% 
% % Calculate the percentage of shoreline change
% shoreline_change_percentage = (erosion_length_pixels / shoreline1) * 100;
% 
% % Display the percentage of shoreline change
% disp(['Percentage of Shoreline Change: ' num2str(shoreline_change_percentage) '%']);


% % Clear the workspace and close all figures
% clear all;
% close all;
% clc;
% 
% % Read the images
% image1 = imread('Miami Beach (1).png'); % Replace with the actual path to your first image
% image2 = imread('Miami Beach (2).png'); % Replace with the actual path to your second image
% 
% % Convert the images to grayscale
% image1_gray = rgb2gray(image1);
% image2_gray = rgb2gray(image2);
% 
% % Apply edge detection using the Canny method
% edges1 = edge(image1_gray, 'Canny');
% edges2 = edge(image2_gray, 'Canny');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Threshold the difference image to obtain the binary mask
% threshold = 0.5; % Adjust the threshold as needed
% binary_mask = difference > threshold;
% 
% % Create a color mask based on the binary mask
% color_mask = cat(3, binary_mask, ~binary_mask, ~binary_mask);
% 
% % Combine the color mask with the original image1
% output = image1;
% output(color_mask) = image2(color_mask);
% 
% % Display the output image
% figure;
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
% imshow(output);
% title('Output Image');
% 
% % Calculate the length of the shoreline in both images
% shoreline1 = sum(edges1(:));
% shoreline2 = sum(edges2(:));
% 
% % Calculate the erosion length in pixels
% erosion_length_pixels = shoreline2 - shoreline1;
% 
% % Display the erosion length in pixels
% disp(['Erosion Length (Pixels): ' num2str(erosion_length_pixels)]);
% 
% % Convert erosion length to the desired scale (e.g., feet)
% % Assuming a conversion factor of 1 pixel = 0.1 feet
% conversion_factor = 0.;
% erosion_length_feet = erosion_length_pixels * conversion_factor;
% 
% % Display the erosion length in feet
% disp(['Erosion Length (Feet): ' num2str(erosion_length_feet)]);
% 
% % Calculate the percentage of shoreline change
% shoreline_change_percentage = (erosion_length_pixels / shoreline1) * 100;
% 
% % Display the percentage of shoreline change
% disp(['Percentage of Shoreline Change: ' num2str(shoreline_change_percentage) '%']);
% 
% % Calculate the average retreat rate in feet per year
% years = 1; % Number of years between the two images
% average_retreat_rate_feet_per_year = erosion_length_feet / years;
% 
% % Display the average retreat rate
% disp(['Average Retreat Rate (Feet/Year): ' num2str(average_retreat_rate_feet_per_year)]);
% 
% % Combine the two images
% combined = imfuse(image1, image2, 'blend', 'Scaling', 'joint');
% 
% % Calculate the absolute difference between the edge images
% difference = abs(double(edges2) - double(edges1));
% 
% % Apply thresholding to identify areas of erosion
% threshold = graythresh(difference);
% erosionMask = imbinarize(difference, threshold);
% 
% % Apply the erosion mask to the combined image
% combinedErosion = combined;
% combinedErosion(repmat(~erosionMask, 1, 1, 3)) = 0;
% 
% % Display the results
% figure;
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
% imshow(combinedErosion);
% title('Combined with Erosion');

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
% subplot(2, 2, 1);
% imshow(image1);
% title('Image 1');
% 
% subplot(2, 2, 2);
% imshow(image2);
% title('Image 2');
% 
% subplot(2, 2, 3);
% imshow(difference);
% title('Difference');
% 
% subplot(2, 2, 4);
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




