type = input('forward-euler or centered-difference? f/c\n', 's');
h = input('step size\n');
interval = input('interval\n');
y = 1;

iters = interval/h;
xvals = [0];
yvals = [y];


for i = 1:iters
    if type == 'f'
        y = y*h*sin(i*h)+y;
    elseif type == 'c'
        if i == 1
	    y = 1;
	else
	    y = 2 * h * y * sin(i * h) + yvals(i-1, 1);
	end
    end
    xvals = [xvals; i*h];
    yvals = [yvals; y];
end

scatter(xvals, yvals);


