# Image analysis frame XCT
Commented scripts to do some image analysis in MATLAB for live radiography.
If you need to do anything fancier, expecially with grains, I highly recommend checking out SPAM: https://www.spam-project.dev/

Most of the code is just cutting smaller domains within the image / a given frame to make it easier to identify the pile boundaries. There is no fancy tracking or anything.

The following commands and approaches are used again and again:
1) Identify boundaries through thresholding, by using regionprops: https://www.mathworks.com/help/images/ref/regionprops.html
A lot of features can be identified if you know a priori their eccentricity, size (approximately), orientation or convexity (solidity).
2) Gaussian blur + Canny edge detection: canny edge detection uses the color gradient to identify edges. For live radiography, a given threshold can identify noise as a boundary. To remove this, Gaussian blur is applied to avoid pixel-sized boundaries.
https://uk.mathworks.com/help/images/ref/imgaussfilt.html
https://uk.mathworks.com/help/images/ref/edge.html
3) Strel & imdilate imclose imfill imerode: dilate boundaries by applying a structuring element (strel), close them depending on the number of pixels required to obtain a close boundary, fill the newly created closed boundary using imfill/imclose.
https://uk.mathworks.com/help/images/ref/strel.html
https://uk.mathworks.com/help/images/ref/imdilate.html
4) If you already know the orientation of a boundary (e.g. pile tip, or plug in the pile being horizontal), you can smooth out the noise by averaging the grayscale value over the pile width and then identifying the boundary as the location where the first derivative of the intensity (gradient) is significantly higher/lower.
