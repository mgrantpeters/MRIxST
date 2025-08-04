% Assume binary mask is already defined
SpD1 = lay128bin == 1;

[imH, imW] = size(SpD1);
outImage = zeros(imH, imW);

% Parameters
numSquares = 10;
sqSize = 8;
minDist = sqSize*4;
halfSq = floor(sqSize / 2);

% Store centers of squares to check distances
squareCenters = [];

% Find valid positions (top-lefts where a full square can fit)
validMask = conv2(double(SpD1), ones(sqSize), 'valid') == sqSize^2;
[yValid, xValid] = find(validMask);

maxItr = 1000;
k = 0;
itr = 0;

while k < numSquares 
    itr = itr + 1;

    idx = randi(length(xValid));
    x1 = xValid(idx);
    y1 = yValid(idx);
    x2 = x1 + sqSize - 1;
    y2 = y1 + sqSize - 1;

    center = [x1 + halfSq, y1 + halfSq];

    % Distance check
    if ~isempty(squareCenters)
        dists = sqrt(sum((squareCenters - center).^2, 2));
        if any(dists < minDist)
            continue
        end
    end

    % Passed all checks
    outImage(y1:y2, x1:x2) = 1;
    squareCenters = [squareCenters; center];
    k = k + 1;
    disp(['Placed square: ' num2str(k)]);
end

% Display results
figure; imshow(lay128bin, []);
title('Original Binary Mask');
figure; imshow(outImage);
title('Generated Squares');
figure(); imshow(double(lay128bin) - outImage,[]);
title('Diff');