
x = linspace(-pi,pi,100);

subplot(1,2,1);
plot(x,sin(x),'r');
axis equal

subplot(1,2,2);
plot(x,sin(x),'b');
axis equal

figure
[ha, pos] = tight_subplot(1,2,[.05 .05],[.05 .05],[.06 .01]);
for ii = 1:2
    axes(ha(ii)); 
    plot(x,sin(ii*x));
    xlabel('x')
    ylabel('y')
end
set(ha(1:2),'XTickLabel',''); set(ha,'YTickLabel','')

