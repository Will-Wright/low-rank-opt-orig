function my_imshow(X)
% my_imshow(X)
%   is similar to the "imshow" command in the image processing
%   toolbox.

imagesc(X)
colormap(gray);
axis image;         % makes is square
axis off;           % gets rid of tick marks
