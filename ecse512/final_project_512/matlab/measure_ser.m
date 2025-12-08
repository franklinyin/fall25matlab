function ser = measure_ser(x_true, x_hard)
% Symbol Error Rate
    ser = mean(x_true ~= x_hard);
end
