A = [3, -4; 1, 2];
y = [-1; 3];
x = [0; 0];

type = input('Use Jacobi iteration or Gauss-Seidel method? j/g:\n', 's');

iters = input('Num of iterations:\n');
disp(A);
disp(y);

sizeA = size(A);
xcur = zeros(sizeA(1,1), 1);
xnxt = xcur;
xPoints = [];
yPoints = [];
iPoints = [];

if type == 'j'
    for i = 1:iters
        for j = 1:sizeA(1,1)
            rowSum = 0;
            for k = 1:sizeA(1, 2)
                rowSum = rowSum + A(j, k) * xcur(k, 1);
	    end
	    rowSum = rowSum - y(j, 1);
            xnxt(j, 1) = xcur(j, 1) - (1/A(j, j)) * rowSum;
        end
	xcur = xnxt;
        xPoints = [xPoints; xcur(1, 1)];
        yPoints = [yPoints; xcur(2, 1)];
        iPoints = [iPoints; i];
    end
end

if type == 'g'
    for i = 1:iters
        %{
        for j = 1:sizeA(1,1)
            rowSum = 0;
            for k = 1:sizeA(1, 2)
                rowSum = rowSum + A(j, k) * xcur(k, 1);
	    end
	    rowSum = rowSum - y(j, 1);
            xcur(j, 1) = xcur(j, 1) - (1/A(j, j)) * rowSum;
        end
	%}
        xcur(1,1) = (-1 + 4 * xcur(2, 1)) / 3;
	xcur(2,1) = (3 - xcur(1, 1)) / 2;
        xPoints = [xPoints; xcur(1, 1)];
        yPoints = [yPoints; xcur(2, 1)];
        iPoints = [iPoints; i];
    end
end

	        

disp(xcur);
scatter3(xPoints, yPoints, iPoints);
