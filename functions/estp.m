function p_est = estp(q,r,y,D)
% activity rate estimator
p_est = 1/D*(norm(y,q)/norm(y,r))^(1/(1/q-1/r));
end