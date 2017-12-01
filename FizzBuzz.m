% solving fizz buzz in matlab
clear print
for i = 1:100
    if i/15 == floor(i/15)
        print{i} = 'FizzBuzz';
    elseif i/3 == floor(i/3)
        print{i} = 'Fizz';
    elseif i/5 == floor(i/5)
        print{i} = 'Buzz';
    else
        print{i} = num2str(i);
    end
end
print

% or
clear print
for i = 1:100
    s = '';
    if mod(i,3) > 0 && mod(i,5) > 0
        s = horzcat(s,num2str(i));
    elseif mod(i,3) == 0
        s = horzcat(s,'Fizz');
    end
    if mod(i,5) == 0
        s = horzcat(s,'Buzz');
    end
    print{i} = s;
end
print

% or
clear print
for i = 1:100
    s = '';
    if xor(mod(i,3)>0,mod(i,5)>0) == 1
        if mod(i,3) == 0
            s = horzcat(s,'Fizz');
        else
            s = horzcat(s,'Buzz');
        end
    else
        s = horzcat(s,num2str(i));
    end
    print{i} = s;
end
print

% or, from Cody
for d= 1:100
    if mod(d,15)==0
        disp('fizz buzz')
    elseif mod(d,3)==0
        disp('fizz')
    elseif mod(d,5)==0
        disp('buzz')
    else
        disp(d);
    end
end