function detSpectrum = draw_spectogram(handles)

warning off;
fftl = handles.FFTLVal;
skip = handles.OverlapVal;

nRecs = round((length(handles.AudioData)/fftl)*(fftl/skip)); %allows for last call to be seen fully

freq = (0:fftl/2)/fftl*handles.SampleFreqVal;
[~,low] = min(abs(freq-handles.StartFreqVal));
[~,high] = min(abs(freq-handles.EndFreqVal));
win = hamming(fftl);

%x=data;
x = [handles.AudioData; zeros(fftl,1)]; 
[x1,x2] = size(x);
if x2>x1
    x=x';
end

detSpectrum = zeros(high-low+1,nRecs);
for j = 1:nRecs
    start = (j-1)*skip+1;
    finish = start+fftl-1;
    q = fft(x(start:finish).*win);
    detSpectrum(:,j) = abs(q(low:high));
end

% If whitening field is on
if(handles.Whiten == 1)
    nFreqBins = high-low+1;
    [spc,~,mu] = whiten3(detSpectrum);
    spc = abs((spc./(mu*ones(1,nRecs))).')';
    
    [spc,~,mu] = whiten3(spc');
    spc=abs(spc'./(ones(nFreqBins,1)*mu'));
    
    cross=ones(3,3);cross(2,2)=4;cross(1,3)=0;
    cross(3,1)=0;cross(1,1)=0;cross(3,3)=0;cross=cross/8;
    spc=conv2(spc,cross,'same');
    detSpectrum=spc;
end
% normalize so brightness can be semi-predictable
detSpectrum = detSpectrum/max(max(detSpectrum));

specLength = size(handles.dim_coords);
if(specLength(2)>1)
    dim = zeros(size(detSpectrum));
    dim(handles.dim_coords(1,1):handles.dim_coords(1,2),:) = 1;
    imagesc(1:size(detSpectrum,2),freq(low:high),(dim.*detSpectrum).^((handles.brightness)*.5));
else
    imagesc(1:size(detSpectrum,2),freq(low:high),(detSpectrum).^((handles.brightness)*.5));
end
set(gca,'CLim',[0,.5+(handles.brightness*.5)])

% Set up axis labels
axis xy;
[specLength,~] = size(detSpectrum);
% set(gca,'XTick',1:floor(nRecs/20):nRecs)
% set(gca,'XTickLabel',0:handles.PlotLengthVal/20:handles.PlotLengthVal)
% 
% set(gca,'YTick',1:20:specLength)
% set(gca,'YTickLabel',handles.StartFreqVal:round(20*(handles.EndFreqVal-handles.StartFreqVal)/specLength):handles.EndFreqVal)
set(gca,'TickDir','out')
if(handles.deliminate_calls==1)
    hold on;
    for iDilim = 1:length(handles.markers)-1
%     while(handles.markers(i)<(handles.PlotLengthVal*handles.SampleFreqVal) && ...
%             i + sum(handles.reverse_vector(1:end-1)) <= length(handles.bt))
        
         plot([(handles.markers(iDilim)-(fftl/2-skip))/skip,...
             (handles.markers(iDilim)-(fftl/2-skip))/skip],...
             [0,handles.EndFreqVal],'w');
        %       h1=datestr(bt(i+sum(reverse_vector(1:end-1)),5),'HH'); h1=str2num(h1);
        %       h2=datestr(bt(i+1+sum(reverse_vector(1:end-1)),4),'HH');h2=str2num(h2);
        %       d1=datestr(bt(i+sum(reverse_vector(1:end-1)),5),'DD'); d1=str2num(d1);
        %       d2=datestr(bt(i+1+sum(reverse_vector(1:end-1)),4),'DD');d2=str2num(d2);
        % if(h2>h1 || d2>d1)
        %         if(h2>h1)
        %         num_hours=h2-h1;
        %         else
        %             num_hours=1;
        %         end
        % plot([(markers(i)-(fftl/2-skip))/skip,(markers(i)-(fftl/2-skip))/skip],[0,sample_freq/2],'r','LineWidth',num_hours+1);
        % end
        % i+sum(reverse_vector(1:end-1))
        
        boxnumber = handles.bt(handles.ViewStart+iDilim-1,3);
        
        boxcolor(1)='r';
        boxcolor(2)='g';
        boxcolor(3)='c';
        boxcolor(4)='b';
        boxcolor(5)='m';
        boxcolor(6)='y';
        boxcolor(7:100)='w';
        
        rectangle('Position',[(handles.markers(iDilim)-(fftl/2-skip))/skip,...
            handles.StartFreqVal,...
            ((handles.markers(iDilim+1)-(fftl/2-skip))/skip)-...
            ((handles.markers(iDilim)-(fftl/2-skip))/skip),...
            round(.01*(handles.EndFreqVal-handles.StartFreqVal))],...
            'FaceColor',boxcolor(boxnumber+1))
        thisExpNum = handles.ViewStart+ iDilim-1;
        if thisExpNum>1 && thisExpNum<=size(handles.bt,1)
            % calculate time since previous detection and round to nearest
            % second.
            dtBT = round(24*60*60*(handles.bt(thisExpNum,4)-handles.bt(thisExpNum-1,4))); 
            text((handles.markers(iDilim)-(fftl/2-skip))/skip,handles.EndFreqVal*.8,...
                sprintf('%ds',dtBT),...
                'FontSize',12,'FontWeight','demi','BackgroundColor','w')
        else
            text((handles.markers(iDilim)-(fftl/2-skip))/skip, handles.EndFreqVal*.8,'unk.',...
                'FontSize',12,'FontWeight','demi','BackgroundColor','w')
        end   
        %       if(bt(i+sum(reverse_vector(1:end-1)),3)==0)
        %           rectangle('Position',[(markers(i)-(fftl/2-skip))/skip,0,...
        % ((markers(i+1)-(fftl/2-skip))/skip)-((markers(i)-(fftl/2-skip))/skip),3],'FaceColor','r')
        %       end
        %       if(bt(i+sum(reverse_vector(1:end-1)),3)==1)
        %           rectangle('Position',[(markers(i)-(fftl/2-skip))/skip,0,((markers(i+1)-(fftl/2-skip))/skip)-((markers(i)-(fftl/2-skip))/skip),3],'FaceColor','g')
        %       end
        %
        %       if(bt(i+sum(reverse_vector(1:end-1)),3) > 1)
        %           rectangle('Position',[(markers(i)-(fftl/2-skip))/skip,0,((markers(i+1)-(fftl/2-skip))/skip)-((markers(i)-(fftl/2-skip))/skip),3],'FaceColor','c')
        %       end
        % iDilim = iDilim+1;
    end
    % rectangle('Position',[markers(i)/skip,0,(markers(i+1)/skip)-(markers(i)/skip),3],'FaceColor','g')
    hold off;
end
ylim([handles.StartFreqVal,handles.EndFreqVal])