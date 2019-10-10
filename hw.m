% compute the hamming weight for every matrix element

function result = hw(input)

    result = arrayfun(@(el) sum(dec2bin(el)=='1'), input);

end