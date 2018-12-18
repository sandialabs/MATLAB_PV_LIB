% pvl_lambertw_test
% Jonathan Allen 2018-12-18

%% set up test
% start with large number of random numbers to find a few that do not
% converge to machine precision
rng('default'); % initialize seed
z = 10 * rand(1e3,1);

% keep track of convergence
w2 = zeros(length(z),36);
converge = true(length(z),36);

%% run pvl_lambertw
% Use asymptotic expansion w = log(z) - log(log(z)) for most z
tmp = log(z + (z == 0));
w = tmp - log(tmp + (tmp == 0));

% Use a series expansion when close to the branch point -1/e
k = (abs(z + 0.3678794411714423216) <= 1.5);
tmp = sqrt(5.43656365691809047*z + 2) - 1; % [1], Eq. 4.22 and text
w(k) = tmp(k);

for k = 1:36
   % Converge with Halley's method ([1], Eq. 5.9), about 5 iterations satisfies
   % the tolerance for most z
   c1 = exp(w);
   c2 = w.*c1 - z;
   w1 = w + (w ~= -1);
   dw = c2./(c1.*w1 - ((w + 2).*c2./(2*w1)));
   w = w - dw;

   w2(:,k) = w;
   converge(:,k) = abs(dw) < 0.7e-16*(2+abs(w));
   
   if all(abs(dw) < 0.7e-16*(2+abs(w)))
%* if all(abs(dw) < eps*(2+abs(w))) %* recommended fix
      break;
   end
end

%% analyze test
figure
imagesc(converge);

% non convergent values
idx = find(~converge(:,35) | ~converge(:,36));
for i = 1:length(idx)
    fprintf('z = %f, w = %f \n',z(idx(i)),w(idx(i)));
end

