% This is one of alpha calculation method: AUT method
% Based on the quad_wright code.
% 2014.7.13
% Yi

function [alpha] = Alpha_AUT(y,lambda,delta_y,delta_lambda)

alpha_hat = 0.999995;

%Decide on suitable affine step-length
duals = [y;lambda];
delta = [delta_y;delta_lambda];
index = delta < 0;
if(any(index))
    alpha_hat = 0.9995*min(-duals(index)./delta(index)); %solves for min ratio (max value of alpha allowed)
end
alpha = min([alpha_hat,0.999995]);
%Check for numerical problems (alpha will go negative)
if(alpha < 0), error('alpha<0'); end

end