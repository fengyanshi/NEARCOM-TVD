clear all
%wave=load('wavecon.n1f');
% time
nhour=120; % 5 days
wave=zeros(nhour,5);
wave(:,1)=[0:nhour-1];
tmp_time=wave(:,1);
wave_time=datenum('November 14, 2008 16:00:00.000 AM')+tmp_time/24.;
wave_date=datevec(wave_time);
for k=1:length(wave_date)
   nyear=int2str(wave_date(k,1));
     if wave_date(k,2)<10
    nmonth=['0' int2str(wave_date(k,2))];
   else
    nmonth=int2str(wave_date(k,2));      
     end
      if wave_date(k,3)<10
    ndate=['0' int2str(wave_date(k,3))];
   else
    ndate=int2str(wave_date(k,3))      ;
      end
     if wave_date(k,4)<10
    nhour=['0' int2str(wave_date(k,4))];
   else
    nhour=int2str(wave_date(k,4))      ;
     end
     if wave_date(k,5)<10
    nminute=['0' int2str(wave_date(k,5))];
   else
    nminute=int2str(wave_date(k,5))      ;
     end
     yeartime(k,:)=[nyear nmonth ndate '.' nhour nminute '00'];
end


%plot(wind_date,wind(:,2));
wave_height(1:length(wave))=0.182;
wave_period(1:length(wave))=2.5;
wave_angle(1:length(wave))=10.0;

for k=1:length(wave)
   yt=yeartime(k,1:15);
   hh=num2str(wave_height(k));
   wp=num2str(wave_period(k));
   wa=num2str(wave_angle(k));
   tot=[yt ' ' hh ' ' wp ' ' wa ' ' '10'];
   wavenew(k,1:length(tot))=tot;
end

fid=fopen('wave_lab.txt','w','n');
fprintf(fid,'%s\n','TPAR')
for k=1:length(wavenew)
fprintf(fid, '%s\n', strtrim(wavenew(k,:)));
end
fclose(fid);

break

%z_model=load('outputz.out');

%   model.z=z_model(:,2:19);
 %  model.time=datenum('June 19, 2005  0:00:00.000 AM')+z_model(:,1)/3600./24.;
%save('model_z_n.mat','model');
