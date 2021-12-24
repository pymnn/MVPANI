function y_out = interp1dec(x,y,x_out)

% Alternative version of interp1q that might be slower, but does not
% require the signal processing toolbox

y_out = zeros(size(x_out));
y_out(x) = y;

vals = diff(x);

for i = 1:length(x)-1
    y_out(x(i):(x(i+1))) = linspace(y(i),y(i+1),vals(i)+1);
end