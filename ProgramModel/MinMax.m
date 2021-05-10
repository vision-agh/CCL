function [min, max] = MinMax(a, b)

if a > b
    max = a;
    min = b;
else
    max = b;
    min = a;
end

end