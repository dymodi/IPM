% This is one of alpha calculation method NTU method
% Based on the NTU Master Paper by Wu Xuepei
% 2014.7.13
% Yi

function [alpha] = Alpha_NTU(y,lambda,delta_y,delta_lambda)

alpha_1 = 1;
alpha_2 = 1;

while(min(y+alpha_1*delta_y)<0)
    alpha_1 = alpha_1 -0.001;
    if(alpha_1<=0)
        break;
    end
end

while(min(lambda+alpha_2*delta_lambda)<0)
    alpha_2 = alpha_2 -0.001;
    if(alpha_2<=0)
        break;
    end
end

alpha = min([max([alpha_1,0.001]),max([alpha_2,0.001])]);

end
