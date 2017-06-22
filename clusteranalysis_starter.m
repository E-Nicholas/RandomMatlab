%cluster analysis script starter

%OK, pair-wise t-tests... using ERPavg (ind,cond,chan,time), final map is
%channel/time for 1 condition pair
ERPavg_copy = ERPavg;
ERPavg_copy(:,17,:,:) = bsxfun(@plus,ERPavg(:,1,:,:),ERPavg(:,2,:,:)); %pure SUM is SUM of pureAud(1) and pureVis(2), Compare to pure AV(3)
ERPavg_copy(:,18,:,:) = bsxfun(@plus,ERPavg(:,6,:,:),ERPavg(:,7,:,:)); %repeat SUM is SUM of A2A(6) and V2V(7), Compare to AV2AV(11)
ERPavg_copy(:,19,:,:) = bsxfun(@plus,ERPavg(:,4,:,:),ERPavg(:,9,:,:)); %switch SUM is SUM of V2A(4) and A2V(9), Compare to A2AV or V2AV (13)


%Classic comparison (AV vs SUM A + V)
cond1 = 15;
cond2 = 17;
resptime = mean(inddata(:,17,1),1);

%Pure comparison (Pure AV vs SUM Pure A + Pure V)
cond1 = 3;
cond2 = 18;
resptime = mean(inddata(:,5,1),1);

%Repeat comparison (Repeat AV vs SUM Repeat A + Repeat V)
cond1 = 11;
cond2 = 19;
resptime = mean(inddata(:,13,1),1);

%Switch comparison (Switch AV vs SUM Switch A + Switch V)
cond1 = 16;
cond2 = 20;
%resptime = mean(inddata(:,18,1),1); 
resptime = 207.6

for i = 1:length(ERPavg_copy(1,1,:,1)) %calculating statistics here
   
   for j = 1:length(ERPavg_copy(1,1,1,:))
       [hvals(i,j),pvals(i,j),civals(i,j,:),allstats{i,j}] = ttest(ERPavg_copy(:,cond1,i,j),ERPavg_copy(:,cond2,i,j));
   end
    
end

for i = 1:length(hvals(:,1)) %making sure there are a series of significant points in a row before reporting any individual blip of significance
    
    for j = 1:length(hvals(1,:))
        if hvals(i,j) == 1
            if j < 6
                hcorrect(i,j) = 0;
            elseif j > length(hvals(1,:))-6
                hcorrect(i,j) = 0;
            elseif hvals(i,j+1) == 1 && hvals(i,j+2) == 1 && hvals(i,j+3) == 1 && hvals(i,j+4) == 1 && hvals(i,j-1) == 1 && hvals(i,j-2) == 1 && hvals(i,j-3) == 1 && hvals(i,j-4) == 1
                hcorrect(i,j) = 1;
            else
                hcorrect(i,j) = 0;
            end
        else
            hcorrect(i,j) = 0;
        end
    end
    
end


%creating plotmap of t-stat values
for i = 1:length(hvals(:,1))
    
   for j = 1:length(hvals(1,:))
      if hvals(i,j) == 1
        plotmap(i,j) = allstats{i,j}.tstat;
      else
          plotmap(i,j) = 0;
      end
   end
    
end


%creating plotmap... I don't know how this is different from above, other
%than it works on hcorrect instead of hvals, it also re-writes the
%variable... I may have left the above loop because I was in a rush and/or
%being sloppy, or maybe it's important, I'm too lazy to test now.
for i = 1:length(hcorrect(:,1))
    
   for j = 1:length(hcorrect(1,:))
      if hcorrect(i,j) == 1
        plotmap(i,j) = allstats{i,j}.tstat;
      else
          plotmap(i,j) = 0;
      end
   end
    
end

%Ok, now reorder the electrodes to something that systematically works its
%way across the head and is good for a linear depiction of channels
neworder = {'H19' 'H18' 'H17' 'H16' 'H15' 'H14' 'H26' 'H27' 'H28' 'H29' 'H30' 'H31' 'H32' 'A9' 'A10' 'A11' 'A12' 'A13' 'A14' 'A15' 'A22' 'A21' 'A20' 'A19' 'A18' 'A17' 'A16' 'A26' 'A27' 'A28' 'A29' 'A30' 'A31' 'A32' 'B6' 'B7' 'B8' 'B9' 'B10' 'B11' 'B12' 'B18' 'B17' 'B16' 'B15' 'B14' 'B13' 'G15' 'G17' 'H7' 'H20' 'G16' 'H5' 'H6' 'H21' 'H24' 'H25' 'G4' 'H4' 'H3' 'H22' 'A7' 'A8' 'F21' 'G3' 'H2' 'H23' 'F20' 'G2' 'H1' 'F1' 'G1' 'A1' 'A2' 'A3' 'A4' 'A5' 'A6' 'A23' 'D1' 'C1' 'D2' 'C2' 'B1' 'C25' 'C3' 'B21' 'B2' 'C24' 'C4' 'B22' 'B3' 'A24' 'A25' 'C23' 'C5' 'B23' 'B20' 'B4' 'B5' 'C22' 'C6' 'B24' 'B19' 'G32' 'H13' 'G22' 'G31' 'H12' 'F30' 'F31' 'G10' 'G11' 'G21' 'G23' 'G30' 'H11' 'F32' 'G9' 'G12' 'G20' 'G24' 'G29' 'H10' 'G8' 'G13' 'G19' 'G25' 'G28' 'H9' 'G7' 'G14' 'G18' 'G26' 'G27' 'H8' 'D6' 'C29' 'C21' 'C7' 'C8' 'B25' 'D7' 'C30' 'C20' 'C15' 'C9' 'B26' 'D8' 'C31' 'C19' 'C16' 'C14' 'C10' 'B27' 'D10' 'D9' 'C32' 'C18' 'C17' 'C13' 'C11' 'B28' 'C12' 'B32' 'B29' 'B31' 'B30' 'G6' 'G5' 'F24' 'F23' 'F22' 'F16' 'F17' 'F18' 'F19' 'F6' 'F5' 'F4' 'F3' 'F2' 'E24' 'E23' 'E22' 'E21' 'E16' 'E17' 'E18' 'E1' 'E20' 'E19' 'E6' 'E7' 'E2' 'E3' 'E4' 'E5' 'D28' 'D3' 'D4' 'D17' 'D18' 'C26' 'D5' 'D16' 'C27' 'C28' 'F29' 'F28' 'F27' 'F26' 'F25' 'F11' 'F12' 'F13' 'F14' 'F15' 'F10' 'F9' 'F8' 'F7' 'E29' 'E30' 'E31' 'E32' 'E28' 'E27' 'E26' 'E25' 'E12' 'E13' 'E14' 'E15' 'E11' 'E10' 'E9' 'E8' 'D32' 'D31' 'D30' 'D29' 'D24' 'D25' 'D26' 'D27' 'D23' 'D22' 'D21' 'D20' 'D19' 'D11' 'D12' 'D13' 'D14' 'D15'};
%given new order above, regions of interest are as follows: start thru G15
%-> Occipital, G15->G32-> Parietal, G32(Temporal Left)->D6->(Temporal
%Right)B30, Central B30->C28, Frontal -> F29
occreg = [1 find(strcmp('G15',neworder))];
parietreg = [find(strcmp('G17',neworder)) find(strcmp('G32',neworder))];
lefttempreg = [find(strcmp('H13',neworder)) find(strcmp('D6',neworder))];
righttempreg = [find(strcmp('C29',neworder)) find(strcmp('B30',neworder))];
centralreg = [find(strcmp('G6',neworder)) find(strcmp('C28',neworder))];
frontalreg = [find(strcmp('F29',neworder)) length(neworder)];

finalorder = [neworder(frontalreg(1):frontalreg(2)) neworder(centralreg(1):centralreg(2)) neworder(lefttempreg(1):lefttempreg(2)) neworder(righttempreg(1):righttempreg(2)) neworder(parietreg(1):parietreg(2)) neworder(occreg(1):occreg(2))];

frontalreg = [1 find(strcmp(neworder(length(neworder)),finalorder))];
centralreg = [frontalreg(2) find(strcmp('C28',finalorder))];
lefttempreg = [centralreg(2) find(strcmp('D6',finalorder))];
righttempreg = [lefttempreg(2) find(strcmp('B30',finalorder))];
parietreg = [righttempreg(2) find(strcmp('G32',finalorder))];
occreg = [parietreg(2) length(finalorder)];

for i = 1:length(chanlocs)
    chans{i} = chanlocs(i).labels;
end
for i = 1:length(finalorder)
   indexorder(i) = find(strcmp(finalorder{i},chans));
end

for i = 1:length(indexorder)
    orderedplot(i,:) = plotmap(indexorder(i),:);
end

plottimes = -100:50:850;
for i = 1:length(plottimes)
    plottpoints(i) = find(abs(bsxfun(@minus,plottimes(i),t))<1);
end
%Time to generate the actual plots
figure
imagesc(orderedplot);
%imagesc(plotmap);
ax = gca;
ax.XTick = plottpoints;
ax.XTickLabel = plottimes;

[cmin cmax] = caxis; %%THIS IS LOGIC TO PRODUCE A CUSTOM COLORMAP...
%from the range of data currently plotted
crange = linspace(cmin,cmax,64);
lowrange = find(crange<-2);
highrange = find(crange>2);

custommap = ones(64,3);
for i = 1:length(lowrange)
    custommap(i,:) = [i/length(lowrange) i/length(lowrange) 1];
end
count = length(highrange);
for i = highrange(1):highrange(end)
    custommap(i,:) = [1 count/length(highrange) count/length(highrange)];
    count = count - 1;
end
%End colormap logic
colormap(custommap) %apply colormap
hold on;
set(gca,'xlim',[26 206]);
ys = get(gca,'ylim');
plot([52 52],[ys(1) ys(2)],'k--')

%spatialdivs = [occreg(1) parietreg(1) lefttempreg(1) righttempreg(1) centralreg(1) frontalreg(1)];
spatialdivs = [frontalreg(1) centralreg(1) lefttempreg(1) righttempreg(1) parietreg(1) occreg(1)];
ax.YTick = spatialdivs;
%ax.YTickLabel = {'Occ' 'Parietal' 'Left Temp' 'Right Temp' 'Central' 'Frontal'};
ax.YTickLabel = {'Frontal' 'Central' 'Left Temp' 'Right Temp' 'Parietal' 'Occipital'};

patch([0 512 512 0],[frontalreg(1) frontalreg(1) frontalreg(2) frontalreg(2)],[0.95 0.95 0.95],'FaceAlpha',0.25,'EdgeColor','none') 
patch([0 512 512 0],[centralreg(1) centralreg(1) centralreg(2) centralreg(2)],[0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none') 
patch([0 512 512 0],[lefttempreg(1) lefttempreg(1) lefttempreg(2) lefttempreg(2)],[0.95 0.95 0.95],'FaceAlpha',0.25,'EdgeColor','none')
patch([0 512 512 0],[righttempreg(1) righttempreg(1) righttempreg(2) righttempreg(2)],[0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none') 
patch([0 512 512 0],[parietreg(1) parietreg(1) parietreg(2) parietreg(2)],[0.95 0.95 0.95],'FaceAlpha',0.25,'EdgeColor','none') 
patch([0 512 512 0],[occreg(1) occreg(1) occreg(2) occreg(2)],[0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none') %occipital

colorbar

resppoint = find(abs(bsxfun(@minus,resptime,t))<1);
plot([resppoint resppoint],ys,'--g');

%plotting fcn

%neworder = {'H19' 'H18' 'H17' 'H16' 'H15' 'H14' 'B18' 'B17' 'B16' 'B15' 'B14' 'B13' 'H26' 'H27' 'H28' 'H29' 'H30' 'H31' 'H32' 'B6' 'B7' 'B8' 'B9' 'B10' 'B11' 'B12' 'A9' 'A10' 'A11' 'A12' 'A13' 'A14' 'A15' 'A26' 'A27' 'A28' 'A29' 'A30' 'A31' 'A32' 'A22' 'A21' 'A20' 'A19' 'A18' 'A17' 'A16' 'G15' 'C22' 'G17' 'C6' 'H7' 'B24' 'H20' 'B19' 'H25' 'B5' 'A8' 'A25' 'A23' 'G16' 'C23' 'H5' 'C5' 'H6' 'B23' 'H21' 'B20' 'H24' 'B4' 'A7' 'A24' 'A6' 'G4' 'C24' 'H4' 'C4' 'H3' 'B22' 'H22' 'B3' 'A5' 'F21' 'C25' 'G3' 'C3' 'H2' 'B21' 'H23' 'B2' 'A4' 'F20' 'D2' 'G2' 'C2' 'H1' 'B1' 'A3' 'F1' 'D1' 'G1' 'C1' 'A2' 'A1' };


