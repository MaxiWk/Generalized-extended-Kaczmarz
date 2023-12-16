function extracted = rand_sample(n,k)
    assert(n>=k, 'Requires n>=k!')
    vec = 1:n;
    extracted = zeros(1,k);
    for i = 1:k
        rand_index = randi(n-i+1,1);
        extracted(i) = vec(rand_index);
        vec(rand_index) = [];
    end
end