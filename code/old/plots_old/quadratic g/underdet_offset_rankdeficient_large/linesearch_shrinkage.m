function [x,z,t]=linesearch_shrinkage(x,z,a,y,c,b)
    % find argmin g(t)=f(z-t*a)+t*y where grad(f) is shrinkage with c
    % we have g'(t)= m*t-b continuous, increasing and piecewise linear, i.e. slope m and intercept -b vary on different intervals
    % find t with 0=g'(t) <=> t=b/m in the appropriate interval
    % if g'(0)<0 (>0) then we must have t>0 (<0)
    % it is assumed that x=shrinkage(z,c)
    
    % start intervall contains initial value t0=0, i.e. initial b=-g'(0)
    if nargin < 6
        b = a'*x-y;
    end
    s = sign(b);
    if s == 0 % <=> g'(0)=0 and we are already done
        t = 0;
        return;
    elseif s == -1
        a = -a; % to ensure that we look for t>0
        b = -b;
    end
    % from now on only entries are needed where a~=0 
    I = (a~=0);
    n = numel(I);    
    z = z(I);
    a = a(I);
    a_square = a.^2; % will be needed for updates of m and b
    csa = c*sign(a);
    % compute kinks
    r = (z+csa)./a; % right kinks
    l = (z-csa)./a; % left kinks
    m = sum(a_square(l>0)) + sum(a_square(r<0)); % initial M
    % keep track of relevant indices
    I = 1:n; 
    % only kinks >= 0 may be crossed
    Ir = I(r>=0);
    Il = I(l>=0);
    r = r(Ir);
    l = l(Il);
    nl = numel(l); % number of left kinks
    % collect all kinks
    kink = [l;r];
    % and sort them
    [kink,ind] = sort(kink); % now the successive intervals are [kink(k),kink(k+1)]
    % keep track of original indices
    I = [Il,Ir];
    I = I(ind);
    N = numel(kink);
    k = 0;
    while (k < N) && (m*kink(k+1) < b) % <=> g'(kink(k+1)) < 0
        k = k+1;
        if ind(k) <= nl % t crosses some left kink
            b = b-kink(k)*a_square(I(k));
            m = m-a_square(I(k));
        else % t crosses some right kink
            b = b+kink(k)*a_square(I(k));
            m = m+a_square(I(k));
        end    
    end
    if (m == 0) && (N > 0)
        t = s*kink(max(k,1));
    else
        t = s*b/m;
    end
    z = z-t*s*a;
    x = min(max(0,z-c),z+c); % x=shrinkage(z,c)
end
