%noise circle stimulus

%image array 1000x1000

im = rand(501,501);
imcopy = im;

ang=0:0.01:2*pi;
r = 40;
xp=r*cos(ang);
yp=r*sin(ang);

xs = round(xp+250);
ys = round(yp+250);

% for i = 1:length(xp)
%     imcopy(xs(i),ys(i)) = 1;
% end

map = zeros(501,501);
for i = 1:length(im(:,1))
    for j = 1:length(im(:,2))
        if ((i-250)^2) + ((j-250)^2) < r^2
            map(i,j) = 1;
        end
    end
end

idxs = find(map==1);
npix = length(idxs);

imcopy = im;
coherence = 0.05;
frac = round(npix*coherence);

randpix = randsample(idxs,frac,false);

imcopy(randpix) = 1;
imshow(imcopy)

figure(1)
coherence = 0.0;
for n = 1:90
      
      im = rand(501,501);
      imcopy = im;

      
      if n > 60
        coherence = 0.0; 
        frac = round(npix*coherence);
        randpix = randsample(idxs,frac,false);
        imcopy(randpix) = 1;
      elseif n > 30
        coherence = 0.5;
        frac = round(npix*coherence);
        randpix = randsample(idxs,frac,false);
        imcopy(randpix) = 1;
      end
      %imwrite(imcopy,strcat('circle_',num2str(n-1),'.bmp'),'bmp')
      imshow(imcopy)
      pause(0.001)
end
