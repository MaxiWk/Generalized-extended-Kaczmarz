x= linspace(-1,1,100);
y=x.^2;  %get something to plot

h1=subplot(1,2,1); %setup subplot1
plot(x,y,'-.'); %plot subplot1
h1_pos = get(h1,'Position'); %get the position data for sublot1.


h2=subplot(1,2,2);  %make subplot2
y = x.^3;
h2_pos=get(h2,'Position'); 

set(h2,'Position',[h2_pos(1) + 0.05, h2_pos(2:end)]) %using position of subplot1 put subplot2next to it.

