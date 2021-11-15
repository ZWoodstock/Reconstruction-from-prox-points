function out = subgrad_proj(x,f,xi)
%Computes the subgradient projection of f at a vector x with level
%set xi. f must output both the evaluation of f and a selection of
%its subgradient, df. f must be convex, proper, and lsc. This
%operator is firmly quasinonexpansive, or "T-class".
[f, df] = f(x);
temp = xi-f;
if temp<0
  out = x + ((temp)/(norm(df(:)).^2))*df;
else
  out = x;
end
end
  
