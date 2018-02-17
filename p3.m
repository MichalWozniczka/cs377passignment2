deltaT = input('Delta t:\n');
deltaX = input('Delta x:\n');

const = deltaT / (deltaX^2);

cur = zeros(10/deltaT, 10/deltaX);

matSize = size(cur);

xCoords = [];
tCoords = [];
mCoords = [];

for i = 1:matSize(1,2)
    cur(1,i) = exp(-4 * (((i*deltaX) - 5)^2));
    xCoords = [xCoords; i * deltaX];
    tCoords = [tCoords; 1 * deltaT];
    mCoords = [mCoords; cur(1, i)];
end

for i = 2:matSize(1,1)
    for j = 2:matSize(1,2)-1
        cur(i,j) = const * (cur(i-1,j-1) + cur(i-1, j+1) - 2*cur(i-1, j)) + cur(i-1, j);
        xCoords = [xCoords; j * deltaX];
        tCoords = [tCoords; i * deltaT];
        mCoords = [mCoords; cur(i, j)];
    end
end

%scatter3(xCoords, tCoords, mCoords);
mesh(cur);
