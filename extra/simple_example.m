x = linspace(-3, 3, 50);
y1 = sin(x);
y2 = exp(x);
y3 = x.^2 + 3*x + 2;

y = [y1; y2; y3];

per_mode = 1;
neurons = 4;

if per_mode == 1
    prediction = zeros(size(y));
    for i=1:size(y,1)
        net = feedforwardnet(neurons, 'trainbr');
        [net, tr] = train(net, x, y(i, :));
        plotperform(tr);
        savelocation = sprintf("figures/simple_case_loss_l1_n%d_per_mode%d", neurons, i);
        saveas(gcf, savelocation, 'epsc')
        disp(['Saved to ', savelocation])
        prediction(i, :) = net(x);
    end
else
    net = feedforwardnet(neurons, 'trainbr');
    [net, tr] = train(net, x, y);
    prediction = net(x);
end

scatter(x, y1)
hold on
scatter(x, y2)
scatter(x, y3)
plot(x, prediction(1, :))
plot(x, prediction(2, :))
plot(x, prediction(3, :))

legend('sin(x)', 'exp(x)', 'x^2 + 3x + 2', 'predicted sin(x)', 'predicted exp(x)', 'predicted x^2 + 3x + 2')
hold off