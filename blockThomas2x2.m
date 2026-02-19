function x = blockThomas2x2(J,b,n)

% Extract 2x2 block diagonals
A = cell(n,1);   % lower blocks
B = cell(n,1);   % diagonal blocks
C = cell(n,1);   % upper blocks
d = cell(n,1);   % rhs blocks

for i = 1:n
    rows = (2*i-1):(2*i);
    B{i} = J(rows,rows);
    d{i} = b(rows);

    if i > 1
        cols = (2*(i-1)-1):(2*(i-1));
        A{i} = J(rows,cols);
    else
        A{i} = zeros(2);
    end

    if i < n
        cols = (2*(i+1)-1):(2*(i+1));
        C{i} = J(rows,cols);
    else
        C{i} = zeros(2);
    end
end

% forward elimination 
for i = 2:n
    M = A{i} / B{i-1};          % no inverse — linear solve
    B{i} = B{i} - M*C{i-1};
    d{i} = d{i} - M*d{i-1};
end

%backward substitution 
x = zeros(2*n,1);
x_block = B{n} \ d{n};

x(2*n-1:2*n) = x_block;

for i = n-1:-1:1
    x_block = B{i} \ (d{i} - C{i}*x(2*i+1:2*i+2));
    x(2*i-1:2*i) = x_block;
end

end