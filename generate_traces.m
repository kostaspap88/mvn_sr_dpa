function [plaintext, traces, m, C, avgC] = generate_traces(bit_size, no_traces, sigma, k_star)

% PRESENT sbox
sbox = [12 5 6 11 9 0 10 13 3 14 15 8 4 7 1 2];

% random plaintext
plaintext = randi(2^bit_size, no_traces, 1) - 1;

% phi(x,k*) = sbox(xor(x,k*))
for i=1:no_traces
    y(i) = sbox(bitxor(plaintext(i), k_star) + 1);
end

%  assume the true leakage function is identity
traces = y' + normrnd(0, sigma, no_traces, 1);

% compute mean and variance 
m = zeros(2^bit_size,1);
C = zeros(2^bit_size,1);
for i=0:2^bit_size-1
    
    current_traces = traces(y==i);
    m(i+1) = mean(current_traces);
    C(i+1) = var(current_traces);

end

avgC = mean(C);


end



