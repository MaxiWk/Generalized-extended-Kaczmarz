function x = sparserand_with_abs_min(n, sp, abs_min_value, abs_max_value, real_setting)

%  Produces an column vector of length n with s non-zero entries at random
%  positions and with normally distributed entries (via randn).

x = zeros(n,1);
p = randperm(n);

%rand_signs = (-1) .^ randi(2,1,sp);
x(p(1:sp)) = abs_min_value + (abs_max_value-abs_min_value) * rand(sp,1);

if ~real_setting
    x(p(1:sp)) = abs_min_value + (abs_max_value-abs_min_value) * (1+1i) * rand(sp,1);
end

end
