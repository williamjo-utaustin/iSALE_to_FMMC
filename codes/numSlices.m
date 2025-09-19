function n = numSlices(slice_size)
    total_angle = 2*pi;
    n = total_angle / slice_size;

    if mod(n,1) ~= 0
        error('Slice size does not divide the circle evenly. Try another size.');
    end
end
