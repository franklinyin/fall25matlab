function [R] = CalculateV_SOR_NonUniform(h, N, StepSize)
% The input StepSize is a 1*4 matrix that specify the mesh scheme, StepSize(i)
% StepSize(i) > 0 && StepSize(i) < 1

a = []; bta = [];

% vector for original mesh scheme
for i = 0:1:0.1/h
    x(i+1) = i;
    y(i+1) = i;
end

% change the scheme around (0.06, 0.04)
for i = 1:1:4
    if StepSize(i) < 0 || StepSize(i) >= 1
        error('Step size should less than 1 and bigger than 0')
    end
end

x(0.06/h) = x(0.06/h) + StepSize(1);
x(0.06/h+2) = x(0.06/h+2) - StepSize(2);
y(0.04/h) = y(0.04/h) + StepSize(3);
y(0.04/h+2) = y(0.04/h+2) - StepSize(4);

for i = 1:1:(0.1/h+1)
    if i == (0.1/h+1)
        a(i) = a(i-1);
        b(i) = b(i-1);
    else
        a(i) = x(i+1) - x(i);
        b(i) = y(i+1) - y(i);
    end
end
a
b
R = 0;
end


[R] = CalculateV_SOR_NonUniform(0.01, 10, [0.5, 0.6, 0.3, 0.4])