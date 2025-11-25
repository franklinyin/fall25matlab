function ser = measure_ser(x_true, x_hard)
%MEASURE_SER Symbol Error Rate between ground truth and detected symbols.
    ser = mean(x_true ~= x_hard);
end
