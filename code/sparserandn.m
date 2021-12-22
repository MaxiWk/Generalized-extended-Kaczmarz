function x = sparserandn(n,sp,real_setting)

%  Produces an column vector of length n with s non-zero entries at random
%  positions and with normally distributed entries (via randn).

x = zeros(n,1);
p = randperm(n);

x(p(1:sp)) = randn(sp,1);

if ~real_setting
    x(p(1:sp)) = x(p(1:sp)) + 1i*randn(sp,1);
end

end
