clc
tic 

year = (1981:2014)';
fid = fopen('Fall_lang_pho_spei_3.txt','w');

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
                
                if mean(eos)<=173
                    continue
                end
                
                tran = [year eos];            
                tran(isnan(tran(:,2)),:) = [];
                tran(tran(:,2)==0,:) = [];
                n = length(tran);
                
                if n < 34   %����Ϊ�˼��㷽�㲻����С��34������ݣ�����39����Ԫ��
                    continue
                end
                
                phen = tran;
                t = [yr squeeze(minT(i,j,:))];

                    %��������theta,��һ��ΪPstart,�ڶ���Ϊ�ۻ���ֵ,������Ϊa�����ĸ�Ϊb
                    params = {
                    {'theta(1)', 13, 7, 18} 
                    {'theta(2)', 10, 0, 300}
                    {'theta(3)', 0.5, 0, 1}
%                     {'theta(3)', 0.05, 0, 0.3}
%                     {'theta(4)', 50, 0, 300}  
%                     {'theta(4)', 0.5, 0, 1} 
                    }; 
                
                    %ģ��
                    data.td = [t squeeze(spei(i,j,:))/3];
                    data.pho = Pho;
                    data.lat = i;
                    data.ydata = phen;
                    [res,chain,sigma2,rm] = mcmcrun(model,data,params,options);

                    %�ҵ�ģ�������С��һ�����
                    [rmse,ml] = min(rm);
                    fpa = chain(ml,:);
                    pred = model.modelfun(data,fpa);
                
                    AICc=2*34*log(rmse)+2*3+2*3*(3+1)/(34-3-1); %��������Ϊ4��
                    r = corr(phen(:,2),pred);
                    fprintf(fid,' %6i %6i   %6.2f %6.2f %6.2f   %6.2f %6.1f %6.1f \n',i,j,fpa,r,rmse,AICc);

        end
    end
end
fclose(fid);
toc  
    