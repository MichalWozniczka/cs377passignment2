A = [3, -4; 1, 2];
y = [-1; 3];

deft = input('Use default matrix? y/n:\n');
type = input('Use Jacobi iteration or Gauss-Seidel method? j/g:\n');

if deft == 'n'
    A = input('Enter matrix A:\n');
    y = input('Enter vector y:\n');
end

iters = input('Num of iterations:\n');
disp(A);
disp(y);

sizeA = size(A);
x = zeros(sizeA(1,1), 1);
xPoints = [];
yPoints = [];
iPoints = [];

if type == 'j'
    for i = 1:iters
        for j = 1:sizeA(1,1)
            rowSum = 0;
            for k = 1:sizeA(1, 2)
                rowSum = rowSum + A(j, k) * x(k, 1);
	    end
	    rowSum = rowSum - y(j, 1);
            x(j, 1) = x(j, 1) - (1/A(j, j)) * rowSum;
        end
        xPoints = [xPoints; x(1, 1)];
        yPoints = [yPoints; x(2, 1)];
        iPoints = [iPoints; i];
    end
end

if type == 'g'
    L = zeros(sizeA(1,1), sizeA(1,2));
    U = zeros(sizeA(1,1), sizeA(1,2));
    for j = 1:sizeA(1,1)
        for k = 1:sizeA(1,2)
	    if j < k
	        U(j, k) = A(j, k);
	    else
	        L(j, k) = A(j, k);
	    end
	end
    end
    Linv = inv(L);
    for i = 1:iters
        x = Linv * (y - U * x);
        xPoints = [xPoints; x(1, 1)];
        yPoints = [yPoints; x(2, 1)];
        iPoints = [iPoints; i];
    end
end
	        

disp(xPoints);
scatter3(xPoints, yPoints, iPoints);
