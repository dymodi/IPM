% This is a script to draw some pictures for evaluation.
% Test differences of Iter under different parameter setting and criterions.
% DingYi. 2014.3.29


% subplot(3,1,1);
% plot(Iter_rec,'LineWidth',2,'Color','r');
% title('The setting of \tau');
% xlabel('\theta = 0.3');
% ylabel('Iterations');
% axis normal;

% subplot(3,1,2);
% plot(Iter_rec,'LineWidth',2,'Color','r');
% xlabel('\theta = 0.5');
% ylabel('Iterations');
% axis normal;

% subplot(3,1,3);
% plot(Iter_rec,'LineWidth',2,'Color','r');
% xlabel('\theta = 0.7');
% ylabel('Iterations');
% axis normal;

% subplot(2,1,1);
% title('Degeneration Method');
% plot(Iter_rec,'LineWidth',2,'Color','r');
% xlabel('\alpha - 0.01 every loop');
% ylabel('Iterations');
% axis normal;

subplot(2,1,2);
plot(Iter_rec,'LineWidth',2,'Color','r');
xlabel('\alpha /2 every loop');
ylabel('Iterations');
axis normal;