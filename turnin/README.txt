Group members:

Michal Wozniczka
mw34748

David Prasse
dap3376


P1:

matlab -nodesktop -nosplash -r "p1"

When prompted, enter j to see the Jacobi iterative solution, or g to see the
	Gauss-Seidel iterative solution.

When prompted, enter the number of iterations to calculate.

The program will output the matrix A, followed by the vector y, and the
	estimated solution to the vector x. A 3D scatterplot will be shown,
	where the x and y values represent the components of the solution
	vector at iteration z.


P3:

matlab -nodesktop -nosplash -r "p3"

When prompted, enter the value for deltaT.

When prompted, enter the value for deltaX.

A 3D mesh will be shown, where the x axis represents distance / deltaX, the y
	axis represents time / delta T, and the z axis represents temperature.


P4:

make
./p4 <name of input file>

The program will convert the given dimacs file into a CSR and its transpose in 
	memory and proceed to calculate both the push- and pull-style pagerank 
	solutions, which will be written into push.txt and pull.txt 
	respectively.
