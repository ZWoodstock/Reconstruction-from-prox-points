%This script performs the algorithm from our numerics section in
%two ways. It runs both versions in tandem on the same problem. For
%general reference, variables appended with "_h" refer to the
%version without affine extrapolation, and variables appended with
%"_k" refer to the version with affine extrapolation.  

%The way that the original solution in 'ex1_data' was created was
%by using a slightly modified version of this script (without any
%references to x_sol). Essentially, we ran the exact same script,
%but terminate when norm(xnew_h - xnew_k) is small.

close all; clear all;
page_output_immediately(1)

load 'ex1_data'
N = length(x_true);
%Constraints for C_1 (affine constraint: x is bandlimited.) To
%determine how many Fourier frequency bins remain, one can call
%nnz(band).
band_proportion = .1;
mask = round(N*band_proportion/2);
band = ((1:N)<=(mask+1) & ((1:N)>=(1)))' + (((1:N)>=N+1-mask) & ((1:N)<=N))';
proj_aff = @(x) bandlim(x,band);
%Reproject for good measure.
x_true = proj_aff(x_true);

F1 = @(x) project_monotone(x);
%Constraints for C_2: Do not take the exact value of the level set
%at x_true; we only take an over-estimate of 1.5x its true value.
[true_tv, ~] = tvl1_loss(x_true);
gamma_2= 1.5*true_tv;
sproj = @(x) subgrad_proj(x,@tvl1_loss, gamma_2);

num_rand_obs = 250;
%The vector E is loaded from 'ex1_data' to ensure reproduction of
%the same results. However, E was created using the following
%commented line, and it can be used to create a new instantiation
%of the problem.
%E = randn(length(x_true),num_rand_obs);

%block_size is dimension of \mathcal{G}_k, i.e., random scalar
%products will be projected onto their
%best-increasing-approximations, 10-at-a-time.
block_size = 10;
num_blocks = num_rand_obs / block_size;
if num_blocks~=round(num_blocks)
  fprintf('ERROR: Inconsistent block structure')
  break
end
rho = zeros(num_rand_obs,1);
current_block = 0;
%In the following preprocessing loop, we normalize each column of E
%and then store the inverse of the squared operator norm (2-norm)
%of each block of E. (for the construction of firmly nonexpansive
%F_k), and compute the projection of our data.
LE_norms2 = zeros(num_blocks,1);
for i = 1:num_blocks
  block_ind = block_size*current_block;
  for j = 1:block_size
      E(:,block_ind+j) = E(:,block_ind+j)/norm(E(:,block_ind+j));
  end
  %Inner product step for rho:
  rho(block_ind+1:block_ind+block_size) =...
    		F1(E(:,block_ind+1:block_ind+block_size)'*x_true);
  %Store the inverse squared norms of each linear block
  LE_norms2(i)=norm(E(:,block_ind+1:block_ind+block_size),2).^(-2);
  current_block = mod(current_block+1,num_blocks);
end

current_block = 0;

%ALGORITHM LOOP

%Initialize variables for loop
min_err = Inf;
min_x  = zeros(size(x_true));
x0 = zeros(size(x_true));

xold_h = x0;
xold_k = x0;
eold = Inf;
err = Inf;
tol = 10^-3;
max_evals = 10000;
max_it = round(max_evals / block_size)
k=0;
r_block_h = zeros(size(x_true));
norm_blocks_h = 0;
r_block_k = zeros(size(x_true));
norm_blocks_k = 0;

%Both algorithms activate 3 constraints at each iteration. The
%affine-extrapolated algorithm only needs omega_{i,n} to be 0.5
%since the affine part is incorporated elsewhere in the algorithm.
%The algorithm without affine extrapolation must have
%omega_{i,n}=1/3 since, it activates the affine constraint as if it
%were any other constraint. 
%omega_{i,n} are constants:
%no affine extrapolation
nn_h = 1./3; 
%affine extrapolation
nn_k = 0.5; 
%Initializing variables for gathering statistics on extrapolation
%parameters (lambda_list_*) and errors (err_*).
%One will see that, for the affine extrapolated algorithm, the
%error decreases faster and the extrapolation parameters are
%larger. These variables will not be of length == iterations. They
%are only gathered every 10 iterations or so (modifiable in the
%loop). The lambda_list_* variables can be used to visualize the
%size of the extrapolation parameter for each algorithm.
lambda_list_h = [];
lambda_list_k = [];
err_h = [];
err_k = [];
err_iters = [];
%terminate when df=norm(xnew_h - xnew_k) gets small or if maximum
%iterations are reached.
df=Inf;

while k<max_it && df>tol
    block_ind = block_size*current_block;
%At each iteration, cycle through each block, one  at a time,
%activating via the firmly nonexpansive operators  F_k:
    r_aff = proj_aff(xold_h);
    r_norm = norm(xold_h-r_aff);
    xp = proj_aff(xold_k);
    %\Id - F_k + p_k:
    %Traditional extrapolated Haguazeau
    r_block_h = xold_h + LE_norms2(current_block+1)*E(:,block_ind+1:block_ind+block_size)*...
            (rho(block_ind+1:block_ind+block_size)...
             - F1(E(:,block_ind+1:block_ind+block_size)'*xold_h));
    norm_blocks_h = norm(xold_h - r_block_h);
    %Affine extrapolated Haugazeau
    r_block_k = xp + LE_norms2(current_block+1)*E(:,block_ind+1:block_ind+block_size)*...
            (rho(block_ind+1:block_ind+block_size)...
             - F1(E(:,block_ind+1:block_ind+block_size)'*xp));
    norm_blocks_k = norm(xp - r_block_k);
    sp_h = sproj(xold_h);
    sn_h = norm(xold_h-sp_h);

    sp_k = sproj(xp);
    sn_k = norm(xp-sp_k);

    d_h = nn_h*(sum(r_block_h,2) + r_aff + sp_h) - xold_h;
    d_k = proj_aff(nn_k*(sum(r_block_k,2) + sp_k)) - xp;

    %Extrapolation step without affine exploit
      if nnz(norm_blocks_h)==0 && r_norm ==0 && sn_h==0
        lambda_h = 1;
      elseif mod(k,3)==2
        lambda_h=0.5*(nn_h*(sn_h.^2+r_norm.^2+sum(norm_blocks_h.^2)))./(norm(d_h).^2);
      else
        lambda_h=(nn_h*(sn_h.^2+r_norm.^2+sum(norm_blocks_h.^2)))./(norm(d_h).^2);
      end     
    xnew_h = Q_haug(x0,xold_h,xold_h + lambda_h*d_h);

    %Extrapolation step with affine exploit.
    if nnz(norm_blocks_k) == 0 && sn_k==0
        xnew_k = Q_haug(x0,xold_k,xp);
    elseif mod(k,3)==2
        lambda_k=0.5*(nn_k*(sn_k.^2+sum(norm_blocks_k.^2)))./(norm(d_k).^2);
        xnew_k = Q_haug(x0,xold_k,xp + lambda_k*d_k);
    else
        lambda_k=(nn_k*(sn_k.^2+sum(norm_blocks_k.^2)))./(norm(d_k).^2);
        xnew_k = Q_haug(x0,xold_k,xp + lambda_k*d_k);
    end
    xold_h = xnew_h;
    xold_k = xnew_k;
    k=k+1;
    current_block = mod(current_block+1,num_blocks);
    if mod(k,100)==0
        %Gather statistics on iteration and print to screen
	err_iters = [err_iters, k];
        %These are the extrapolation parameters for both algorithms.
        %See for yourself that the extrapolation parameter is larger for the 
        %affine extrapolated algorithm (*_k), over the standard extrapolated
        %algorithm (*_h).
        lambda_list_h = [lambda_list_h, lambda_h];
        lambda_list_k = [lambda_list_k, lambda_k];
        %Errors
        err_h = [err_h, norm(xnew_h - x_sol)];
        err_k = [err_k, norm(xnew_k - x_sol)];
        %Difference between iterates. df eventually goes to zero, since both
        %versions of the algorithm are solving the same problem.
        df = norm(xnew_h - xnew_k);
        fprintf('%2.1f %%\n',100*k/max_it)
    end
end

figure(1)
subplot(3,1,1)
plot(x_sol,'k')
axis([1 N])
title('Solution')
subplot(3,1,2)
plot(xnew_k,'b')
axis([1 N])
title('Recovery (with affine extrap.)')
subplot(3,1,3)
plot(xnew_h,'color',[0,104,67]./255)
axis([1 N])
title('Recovery (without affine extrap.)')
