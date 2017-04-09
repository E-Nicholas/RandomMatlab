function fibonacci(seed1,seed2,nthplace)

%cd('C:\Users\Eric\Desktop\Everything\LifeDocs')

%computing linear recurrance relation sequence to n'th place using form
%F(n) = F(n-1) + F(n-2)

n = nthplace - 2;

%%
sequence = [seed1,seed2];
for i = 1:n
    n2 = sequence(i);
    n1 = sequence(i+1);
    sequence(i+2) = n2 + n1;
end
sequence

lastarea = sequence(end)*sequence(end);
sumprev = sum(sequence(1:end-1).*sequence(1:end-1));

ratio = lastarea/sumprev
% figure; hold on;
% for i = 1:length(sequence)
%     pos = [0 0 sequence(length(sequence)+1-i) sequence(length(sequence)+1-i)];
%     rectangle('Position',pos,'Facecolor',rand(1,3));
% end
% axis equal

%approx golden spiral
figure; hold on;

for i = 1:length(sequence)
    if i == 1
        pos = [0 0 sequence(i)  sequence(i)];
        rectangle('Position',pos,'Facecolor',rand(1,3));
        
        t = linspace(pi, pi/2, 100);
        r = sequence(i);
        x = r * cos(t) + (pos(1) + sequence(i));
        y = r * sin(t) + pos(2);
        plot(x,y,'k','LineWidth',2)
        
        %text(pos(1)+(sequence(i)/2),pos(2)+(sequence(i)/2),num2str(sequence(i)),'FontSize',10);
        
        plotto = 'down';
    else
        switch plotto
            case 'down'
                pos(1) = pos(1);
                pos(2) = pos(2) - sequence(i);
                pos(3) = sequence(i);
                pos(4) = sequence(i);
                rectangle('Position',pos,'Facecolor',horzcat(1,rand(1,2)));
                
                t = linspace(-pi/2, -pi, 100);
                r = sequence(i);
                x = r * cos(t) + (pos(1) + sequence(i));
                y = r * sin(t) + (pos(2) + sequence(i));
                plot(x,y,'k','LineWidth',2)
                
                %text(pos(1)+(sequence(i)/3),pos(2)+(sequence(i)/2),num2str(sequence(i)),'FontSize',15);
                
                plotto = 'right';
            case 'right'
                pos(1) = pos(1) + sequence(i-1);
                pos(2) = pos(2);
                pos(3) = sequence(i);
                pos(4) = sequence(i);
                rectangle('Position',pos,'Facecolor',horzcat(rand(1,2),1));
                
                t = linspace(3/2*pi, 2*pi, 100);
                r = sequence(i);
                x = r * cos(t) + pos(1);
                y = r * sin(t) + (pos(2) + sequence(i));
                plot(x,y,'k','LineWidth',2)
                
                %text(pos(1)+(sequence(i)/3),pos(2)+(sequence(i)/2),num2str(sequence(i)),'FontSize',15);
                
                plotto = 'up';
            case 'up'
                pos(1) = pos(1) - sequence(i-2);
                pos(2) = pos(2) + sequence(i-1);
                pos(3) = sequence(i);
                pos(4) = sequence(i);
                rectangle('Position',pos,'Facecolor',horzcat(rand,1,rand));
                
                t = linspace(pi/2, 0, 100);
                r = sequence(i);
                x = r * cos(t) + pos(1);
                y = r * sin(t) + pos(2);
                plot(x,y,'k','LineWidth',2)
                
                %text(pos(1)+(sequence(i)/3),pos(2)+(sequence(i)/2),num2str(sequence(i)),'FontSize',15);
                
                plotto = 'left';
            case 'left'
                pos(1) = pos(1) - sequence(i);
                pos(2) = pos(2) - sequence(i-2);
                pos(3) = sequence(i);
                pos(4) = sequence(i);
                rectangle('Position',pos,'Facecolor',rand(1,3));
                
                t = linspace(pi, pi/2, 100);
                r = sequence(i);
                x = r * cos(t) + (pos(1) + sequence(i));
                y = r * sin(t) + pos(2);
                plot(x,y,'k','LineWidth',2)
                
                %text(pos(1)+(sequence(i)/3),pos(2)+(sequence(i)/2),num2str(sequence(i)),'FontSize',15);
                
                plotto = 'down';
        end
    end
end
axis equal
end