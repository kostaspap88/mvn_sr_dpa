function M = prediction_function(x, k)

global estimation_choice

% PRESENT sbox
sbox = [12 5 6 11 9 0 10 13 3 14 15 8 4 7 1 2];

% ID case
if (estimation_choice == 1) || (estimation_choice == 2)
    M = sbox(bitxor(x,k) + 1);
end

% HW case
if estimation_choice == 3
   M = hw(sbox(bitxor(x,k) + 1));
end

end