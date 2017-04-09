% taking a string character by character and converting them to...
% morse code sound output. I would like to expand on the function so it can
% output that sound as a wav (pretty easy), as well as other stuff, I dont
% know...

function morsetranslate
alphanumlut = ' ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.,?';
morselut = {{' '} {'.-' 0 1} {'-...' 1 0 0 0} {'-.-.' 1 0 1 0} {'-..' 1 0 0} {'.' 0} {'..-.' 0 0 1 0},...
    {'--.' 1 1 0} {'....' 0 0 0 0} {'..' 0 0} {'.---' 0 1 1 1} {'-.-' 1 0 1} ['.-..' 0 1 0 0],...
    {'--' 1 1} {'-.' 1 0} {'---' 1 1 1} {'.--.' 0 1 1 0} {'--.-' 1 1 0 1} {'.-.' 0 1 0} {'...' 0 0 0},...
    {'-' 1} {'..-' 0 0 1} {'...-' 0 0 0 1} {'.--' 0 1 1} {'-..-' 1 0 0 1} {'-.--' 1 0 1 1} {'--..' 1 1 0 0},...
    {'.-' 0 1} {'-...' 1 0 0 0} {'-.-.' 1 0 1 0} {'-..' 1 0 0} {'.' 0} {'..-.' 0 0 1 0} {'--.' 1 1 0},...
    {'....' 0 0 0 0} {'..' 0 0} {'.---' 0 1 1 1} {'-.-' 1 0 1} {'.-..' 0 1 0 0} {'--' 1 1} {'-.' 1 0},...
    {'---' 1 1 1} {'.--.' 0 1 1 0} {'--.-' 1 1 0 1} {'.-.' 0 1 0} {'...' 0 0 0} {'-' 1} {'..-' 0 0 1},...
    {'...-' 0 0 0 1} {'.--' 0 1 1} {'-..-' 1 0 0 1} {'-.--' 1 0 1 1} {'--..' 1 1 0 0} {'-----' 1 1 1 1 1},...
    {'.----' 0 1 1 1 1} {'..---' 0 0 1 1 1} {'...--' 0 0 0 1 1} {'....-' 0 0 0 0 1} {'.....' 0 0 0 0 0},...
    {'-....' 1 0 0 0 0} {'--...' 1 1 0 0 0} {'---..' 1 1 1 0 0} {'----.' 1 1 1 1 0} {'.-.-.-' 0 1 0 1 0 1},...
    {'--..--' 1 1 0 0 1 1} {'..--..' 0 0 1 1 0 0}};

%generating sound output
Fs = 44000;
timeVec = (0:1/Fs:0.05);
ltimeVec = (0:1/Fs:0.15);
vol = 1;
charspace = zeros(length(timeVec),1)';
wordspace = zeros((length(timeVec)*7),1)';
letspace  = zeros((length(timeVec)*3),1)';
wave = sin(2 * pi * timeVec * 700);
longwave = sin(2 * pi * ltimeVec * 700);

scrz = get( 0, 'Screensize' );
guisize = [scrz(3)/2 scrz(4)/2 scrz(3)/4 scrz(4)/4];
guifig = figure('Toolbar','none','Menubar','none','Position',guisize,'NumberTitle',...
    'off','Name','Morse Maker');
textbox = uicontrol(guifig,'Style','edit',...
    'Position',[5 guisize(4)/5 guisize(3)-10 guisize(4)/1.25],...
    'Max',2,'HorizontalAlignment','left','Callback',@thetext);
transbutt = uicontrol(guifig,'Style','pushbutton','String','Play Morse','Position',...
    [guisize(3)/3 5 guisize(3)/4 guisize(4)/6],'Callback',@playmorse);

    function thetext(hobj,~)
        inputstring = get(hobj,'String');
    end

    function playmorse(hobj,~)
        output = [];
        outwave = [];
        inputstring = get(textbox,'String');
        for i = 1:length(inputstring)
            char = inputstring(i);
            idx = find(alphanumlut==char);
            if idx == 1
                output = sprintf('%s   ',output);
                outwave = [outwave, wordspace];
            else
                output = sprintf('%s %s',output,morselut{idx}{1});
                for j = 2:length(morselut{idx})
                    if j == 2
                        outwave = [outwave, letspace];
                    end
                    if morselut{idx}{j} == 0
                        outwave = [outwave, wave, charspace];
                    elseif morselut{idx}{j} == 1
                        outwave = [outwave, longwave, charspace];
                    end
                    
                end
            end
        end
        sound(outwave,Fs)
        figure;
        plot(outwave)
    end

end