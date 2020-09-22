function tau = occurence_ratio(plaintext, bit_size)

tau = histcounts(plaintext, 0:2^bit_size) ./ size(plaintext,1);

end