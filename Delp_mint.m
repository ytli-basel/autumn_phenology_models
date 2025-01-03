clc
tic 
year = (1981:2014)';
fid = fopen('Delp(1-f(P))-mint.txt','w');

%%%%%����ģ�����ݡ��������ṹ
%����ģ�⺯��
model.modelfun = @(data,theta) predF(data,theta);
%��������
model.ssfun = @(theta,data) rmseF(theta,data);
%����ģ��Ѱ�Ŵ���
options.nsimu = 5000;
%���á�ģ�������ʾ������0��ʾ����ʾ
options.waitbar = 0;
%�������в���ʾ�����ӡ��Ϣ
options.verbosity = 0;
%ģ�⣬resΪģ��������Ϣ�Ľṹ������chainΪģ�������rmΪģ�����

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
                
                if n < 34   %����Ϊ�˼��㷽�㲻����С��34������ݣ�����39����Ԫ��
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
                
                %��������theta,��һ��Ϊt0,�ڶ���Ϊ������ֵ,������Ϊ��ˮ������ֵ�����ĸ�Ϊ��ˮ�ۻ������������Ϊ�ۻ���ʼʱ��
                params = {
                {'theta(1)', (ltb+utb)/2, ltb, utb} 
                {'theta(2)', (maxdl+mindl)/2.0, mindl, maxdl}
                {'theta(3)', 1, 0, 2}
                {'theta(4)', 1, 0, 2}
                {'theta(5)', 3000, 0, 10000}        
                }; 
            
                %ģ��
                data.tem = t;
                data.pho = po;
                data.ydata = phen;
                [res,chain,sigma2,rm] = mcmcrun(model,data,params,options);

                %�ҵ�ģ�������С��һ�����
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
    