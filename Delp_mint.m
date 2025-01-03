clc
tic 
year = (1981:2014)';
fid = fopen('Delp(1-f(P))-mint.txt','w');

%%%%%设置模型数据、参数、结构
%定义模拟函数
model.modelfun = @(data,theta) predF(data,theta);
%定义误差函数
model.ssfun = @(theta,data) rmseF(theta,data);
%定义模拟寻优次数
options.nsimu = 5000;
%设置“模拟进度显示棒”，0表示不显示
options.waitbar = 0;
%让命令行不显示输出打印信息
options.verbosity = 0;
%模拟，res为模拟结果和信息的结构变量，chain为模拟参数，rm为模拟误差

for i = 1:60
    for j = 1:720
        if nsdmask(i,j) == 1
                disp([('Processing the '),strcat(num2str(i),' ',num2str(j)),(' file')]);
                eos = squeeze(phend(i,j,:));
                eos = round(eos);
                tran = [year eos];            
                tran(isnan(tran(:,2)),:) = [];
                tran(tran(:,2)==0,:) = [];
                n = length(tran);
                
                if n < 34   %这里为了计算方便不考虑小于34年的数据（仅有39个像元）
                    continue
                end
                
                phen = tran;
                t = [yr squeeze(minT(i,j,:))];
                
                tem7 = zeros(34,1);
                temeos = zeros(34,1);
                for ii = 1981:2014
                    inter = t(t(:,1)==ii,2);
                    tem7(ii-1980) = mean(inter(182:213));
                    temeos(ii-1980) = inter(round(mean(phen(:,2))));
                end
                ltb = mean(temeos);
                utb = mean(tem7);
                
                peak = round(Ave_peak(i,j));
                po = Pho(:,i);
                maxdl = max(po);
                mindl = po(min(phen(:,2)));
                
                %参数设置theta,第一个为t0,第二个为积温阈值,第三个为降水总量阈值，第四个为降水累积天数，第五个为累积起始时间
                params = {
                {'theta(1)', (ltb+utb)/2, ltb, utb} 
                {'theta(2)', (maxdl+mindl)/2.0, mindl, maxdl}
                {'theta(3)', 1, 0, 2}
                {'theta(4)', 1, 0, 2}
                {'theta(5)', 3000, 0, 10000}        
                }; 
            
                %模拟
                data.tem = t;
                data.pho = po;
                data.ydata = phen;
                [res,chain,sigma2,rm] = mcmcrun(model,data,params,options);

                %找到模拟误差最小的一组参数
                [rmse,ml] = min(rm);
                fpa = chain(ml,:);
                pred = model.modelfun(data,fpa);
     
                AICc=2*n*log(rmse)+2*5+2*5*(5+1)/(34-5-1);
                r = corr(phen(:,2),pred);
                fprintf(fid,' %6i %6i   %6.2f %6.2f %6.2f %6.2f %6.2f   %6.2f %6.1f %6.1f \n',i,j,fpa,r,rmse,AICc);

        end
    end
end
fclose(fid);
toc  
    