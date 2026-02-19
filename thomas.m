function x = thomas(a, b, c, d)
% a = subdiagonal (length n-1)
% b = diagonal     (length n)
% c = superdiag    (length n-1)
% d = RHS          (length n)

n = length(b);

% Forward elimination
for i = 2:n
    w = a(i-1) / b(i-1);
    b(i) = b(i) - w * c(i-1);
    d(i) = d(i) - w * d(i-1);
end

% Back substitution
x = zeros(n,1);
x(n) = d(n) / b(n);

for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end

end