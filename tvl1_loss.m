function [f, df] = tvl1_loss(x);
  %Outputs the function evaluation and subgradient selection of
  %     \|TV(x)\|_{1} 
N = length(x);
df = zeros(N,1);
df(1:(N-1)) = x(2:N) - x(1:(N-1));
f = sum(abs(df));
%subgradient of the componentwise 1-norm is just the sign of the
%input
df = sign(df);
NNZ = nnz(df);
if NNZ~=0
%apply adjoint for df:
    df(N) = df(N-1);
    df(1)=-df(1);
    df(2:(N-1)) = df(1:(N-2)) - df(2:(N-1));
else %otherwise set subdifferential to zero.
    df = zeros(N,1);
end
end
