function [H, g, r] = GLR(r_new, r)


mu_0 = 0;

r = [r(2:end) r_new];


%Using equation 7.35 from the book

%calculating g from a running window.
g = 1/(2*sigma^2*M)*max(sum(r-mu_0))^2;

if (g > h)
    H = 1;
else
    H = 0;
end







