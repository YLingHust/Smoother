function MaxVol=cal_MaxVol()

    datain = csvread('wind_power_E101-44912-2011-05-5min.csv',1,0);
	
    power_h=reshape(datain(:,5),12,length(datain)/12);
	

    var_powerh=var(power_h,1,1);
	
    csvwrite('var_powerh.csv',var_powerh');
	
    figure;
	
    plot(var_powerh);
	
    xlabel('Hour');
	
    ylabel('Variance');
	
    title(['Variance of the power','(Month)'],'Fontsize',16);
	
    saveas(gcf,'Var_power_month','jpg');

    vmin=min(var_powerh);
	
    vmax=max(var_powerh);
	
    vx_month=linspace(vmin,vmax,100);
	
    f=ksdensity(var_powerh);
	
    figure;
	
    plot(vx_month,f,'Linewidth',2);
	
    xlabel('Hour');
	
    ylabel('Variance');

    [fd,x]=ecdf(var_powerh);
	
    figure;
	
    plot(x,fd,'Linewidth',2);
	
    xlabel('Hour');
	
    ylabel('Variance');
	
    [err,loc]=min(abs(fd-0.95));
	
    thres=x(loc);
	
    devia=zeros(1,length(var_powerh));
	
    h_reg=zeros(1,length(var_powerh));
	
    h_count=1;
	
    for ii=1:length(var_powerh)
	
        if(var_powerh(ii)<thres)
		
            max_devia=abs(max(power_h(:,ii))-mean(power_h(:,ii)));
			
            min_devia=abs(min(power_h(:,ii))-mean(power_h(:,ii)));
			
            devia(ii)=max(max_devia,min_devia);
			
            h_reg(h_count)=ii;
			
            h_count=h_count+1;
			
        end
		
    end
    devia(h_count:end)=[];
	
    h_reg(h_count:end)=[];
	
    figure;
	
    plot(devia);
	
    [fdevia,xdevia]=ecdf(devia);
	
    figure;
	
    plot(xdevia,fdevia);
	
    [err,loc1]=min(abs(fdevia-0.99));
	
    MaxVol=xdevia(loc1);

    S=zeros(size(power_h));
	
    Ave_powerH=zeros(size(power_h));
	
    for is=1:length(var_powerh)
	
        if (any(h_reg==is))
		
            [S(:,is), errs]=findS(power_h(:,is),MaxVol);
			
        else
		
            S(:,is)=zeros(12,1);
			
        end
		
        Ave_powerH(:,is)=power_h(:,is)+S(:,is);
		
    end
	
    Ave_power=reshape(Ave_powerH,1,length(datain));
	
    csvwrite('greenPower-05-91011_X10-smooth.csv',[datain(:,1),Ave_power']);
	
    Ave_num=ceil(Ave_power/300);
	
    Ave_num_ch=[abs(Ave_num(2:end)-Ave_num(1:end-1)),0];
	
    Ave_num_ch_day=sum(reshape(Ave_num_ch,length(Ave_num_ch)/3,3));


%   load('data_reg.mat');

%   load1=csvread('Power_RequestRate_clarknet_access_log_Aug28_2.csv',1,0);

    load1=csvread('cluster_datacenter_power_2011_05_five_minute.csv',1,0);
	
    load1(:,3)=load1(:,3)/1.5-300*0.05*11000;
	
%   load2=csvread('Power_RequestRate_UofS_access_log_1995_6_1-7.csv',1,0);

%   load2(:,3)=load2(:,3)/1.5-300*0.05*11000;

%   load3=csvread('Power_RequestRate_NASA_access_log_Jul95_1995_07_01-07.csv',1,0);

%   load3(:,3)=load3(:,3)/1.5-300*0.05*11000;

%   load4=csvread('Power_RequestRate_calgary_access_log_1994_10_25-31.csv',1,0);

%   load4(:,3)=load4(:,3)/1.5-300*0.05*11000;

%   load5=csvread('Power_RequestRate_UCB-home-IP-846890339-847313219_log.csv',1,0);

%   load5(:,3)=load5(:,3)/1.5-300*0.05*11000;


    x1=zeros(3,length(load1));
	
    x10=[0;0;0];
	
    x10_ave=[0;0;0];
	
    SW=0;
	
    SW_ave=0;
	
    for i1=1:length(load1)
	
        int_Pwind=interp1((0:length(datain)-1)*5,datain(:,5),load1(:,1),'linear');
		
        int_Pwind_ave=interp1((0:length(datain)-1)*5,Ave_power,load1(:,1),'linear');
		
        A1=[300,300,150;300,0,0;1,1,1];
		
        b1=[load1(i1,3),int_Pwind(i1),11000]';
		
        b1_ave=[load1(i1,3),int_Pwind_ave(i1),11000]';
		
        x1=A1\b1;
		
        x1_ave=A1\b1_ave;
		
        if((x1(1)-x10(1))*(x1(2)-x10(2))<0)
		
            SW=SW+min(abs(x1(1)-x10(1)),abs(x1(2)-x10(2)));
        end
		
        x10=x1;
		
        if((x1_ave(1)-x10_ave(1))*(x1_ave(2)-x10_ave(2))<0)
		
            SW_ave=SW_ave+min(abs(x1_ave(1)-x10_ave(1)),abs(x1_ave(2)-x10_ave(2)));
			
        end
		
        x10_ave=x1_ave;
		
    end
	
    fprintf('SW=%d, SW_ave=%d\n',SW,SW_ave);
	
    fprintf('MaxVol=%d\n',MaxVol);


end


