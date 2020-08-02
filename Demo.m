% This code test NST+HT+f-FB vary s/M
clear all;
close all;
clc;
warning 'off'
addpath('NST+HT+f-FB');
N =1000;      % %  Size of the vector in the 'original' input space
M =500;        %      Size of the measurement vector
MaxNbIter   = M;
tolSuccess  = 1e-4; % Tolerance on the relative error for a recovered vector to be consider a success y-Ax

rho=1:10:281;
% rho=1:5:251;
Mn=length(rho);

Trials=500;

Success_s=zeros(Mn,1);
Success_6x=zeros(Mn,1);
Success_10x=zeros(Mn,1);
Success_14x=zeros(Mn,1);
Success_x2=zeros(Mn,1);


Time_s=zeros(Mn,1);
Time_6x=zeros(Mn,1);
Time_10x=zeros(Mn,1);
Time_14x=zeros(Mn,1);
Time_x2=zeros(Mn,1);


NbIter_s=zeros(Mn,1);
NbIter_6x=zeros(Mn,1);
NbIter_10x=zeros(Mn,1);
NbIter_14x=zeros(Mn,1);
NbIter_x2=zeros(Mn,1);


h=waitbar(0,'please wait');
for i=1:Mn
    for j=1:Trials
        A = randn(M,N); % generate random Gaussian matrix
        A = A./repmat(sqrt((sum(A.^2,1))),[M,1]);  % normalized
        s =rho(i);    % sparsity value
        q = randperm(N);
        x = zeros(N,1);
        x(q(1:s))=randn(s,1);   % Gaussian random signal
        y = A*x;  % Generate the observations
        
        y0=A'*((A*A')\y);
        phi=A'*((A*A') \A);  %It is calculated off-line for NST
        %-----------------------------------------------------------------------------
        % NST+HT+f-FB f(x)=s
        tic;
        [xstar1,NbIter1] =NstHtFb(y,y0,A,phi,s,MaxNbIter,tolSuccess);
        Time1 = toc;
        MseNstHtFb_s = norm(x - xstar1)/norm(x);
        if MseNstHtFb_s <=1e-4
            Success_s(i)=Success_s(i)+1;
        end
        Time_s(i)  =Time_s(i)+Time1;
        NbIter_s(i)=NbIter_s(i)+NbIter1;
        %-----------------------------------------------------------------------------
        
        %-----------------------------------------------------------------------------
        % NST+HT+f-FB f(x)=6x;
        tic;
        [xstar2,NbIter2] =NstHtFb_6x(y,y0,A,phi,MaxNbIter,tolSuccess);
        Time2 = toc;
        MseNstHtFb_6x = norm(x - xstar2)/norm(x);
        if MseNstHtFb_6x <=1e-4
            Success_6x(i)=Success_6x(i)+1;
        end
        Time_6x(i)=Time_6x(i)+Time2;
        NbIter_6x(i)=NbIter_6x(i)+NbIter2;
        %-----------------------------------------------------------------------------
        
        %-----------------------------------------------------------------------------
        % NST+HT+f-FB f(x)=10x;
        tic;
        [xstar3,NbIter3] =NstHtFb_10x(y,y0,A,phi,MaxNbIter,tolSuccess);
        Time3 = toc;
        MseNstHtFb_10x = norm(x - xstar3)/norm(x);
        if MseNstHtFb_10x <=1e-4
            Success_10x(i)=Success_10x(i)+1;
        end
        Time_10x(i)=Time_10x(i)+Time3;
         NbIter_10x(i)=NbIter_10x(i)+NbIter3;
        %-----------------------------------------------------------------------------
        
        %-----------------------------------------------------------------------------
        % NST+HT+f-FB f(x)=14x;
        tic;
        [xstar4,NbIter4] =NstHtFb_14x(y,y0,A,phi,MaxNbIter,tolSuccess);
        Time4 = toc;
        MseNstHtFb_14x = norm(x - xstar4)/norm(x);
        if MseNstHtFb_14x <=1e-4
            Success_14x(i)=Success_14x(i)+1;
        end
        Time_14x(i)=Time_14x(i)+Time4;
         NbIter_14x(i)=NbIter_14x(i)+NbIter4;
        %-----------------------------------------------------------------------------
        
        %-----------------------------------------------------------------------------
        % NST+HT+f-FB f(x)=x^2;
        tic;
        [xstar5,NbIter5] =NstHtFb_x2(y,y0,A,phi,MaxNbIter,tolSuccess);
        Time5 = toc;
         MseNstHtFb_x2 = norm(x - xstar5)/norm(x);
        if MseNstHtFb_x2 <=1e-4
            Success_x2(i)=Success_x2(i)+1;
        end
        Time_x2(i)=Time_x2(i)+Time5;
         NbIter_x2(i)=NbIter_x2(i)+NbIter5;
        %-----------------------------------------------------------------------------
        
    end
    waitbar(i/Mn,h);
end

NewFolderName = ['fig1'   '_' datestr(now,'mm-dd-yyyy-HH-MM')];
mkdir(NewFolderName)
SourceFileName = [pwd '\' mfilename '.m'];
DestinationFileName = [pwd '\' NewFolderName '\' mfilename '-' datestr(now,'mm-dd-yyyy-HH-MM') '.m.backup'];
copyfile(SourceFileName,DestinationFileName);

cd(NewFolderName)

save('Success_s.mat','Success_s');
save('Success_6x.mat','Success_6x');
save('Success_10x.mat','Success_10x');
save('Success_14x','Success_14x');
save('Success_x2.mat','Success_x2');


save('Time_s.mat','Time_s');
save('Time_6x.mat','Time_6x');
save('Time_10x.mat','Time_10x');
save('Time_14x.mat','Time_14x');
save('Time_x2.mat','Time_x2');


save('NbIter_s.mat','NbIter_s');
save('NbIter_6x.mat','NbIter_6x');
save('NbIter_10x.mat','NbIter_10x');
save('NbIter_14x.mat','NbIter_14x');
save('NbIter_x2.mat','NbIter_x2');

figure;
plot(rho, Success_s/Trials,'-kp','LineWidth',2)
hold on;
plot(rho, Success_6x/Trials,'-m+','LineWidth',2)
plot(rho, Success_10x/Trials,'-b*','LineWidth',2)
plot(rho, Success_14x/Trials,'-gs','LineWidth',2)
plot(rho, Success_x2/Trials,'-rd','LineWidth',2)
xlabel('Sparsity','fontsize',11,'FontWeight','bold');
ylabel('Frequency of exact recovery','fontsize',11,'FontWeight','bold');
legend('f(k)=s','f(k)=6k','f(k)=10k','f(k)=14k','f(k)=k^2');
set(gca,'fontsize',11,'FontWeight','bold');
title('Gaussian sparse vectors');
BMPName = 'Frequency of exact recovery';
print(gcf,'-dbmp',BMPName)


figure;
plot(rho, Time_s/Trials,'-kp','LineWidth',2)
hold on;
plot(rho, Time_6x/Trials,'-m+','LineWidth',2)
plot(rho, Time_10x/Trials,'-b*','LineWidth',2)
plot(rho, Time_14x/Trials,'-gs','LineWidth',2)
plot(rho, Time_x2/Trials,'-rd','LineWidth',2)
xlabel('Sparsity','fontsize',11,'FontWeight','bold');
ylabel('Running time (s)','fontsize',11,'FontWeight','bold');
legend('f(k)=s','f(k)=6k','f(k)=10k','f(k)=14k','f(k)=k^2');
set(gca,'fontsize',11,'FontWeight','bold');
title('Gaussian sparse vectors');
BMPName = 'Running time';
print(gcf,'-dbmp',BMPName)


figure;
plot(rho, NbIter_s/Trials,'-kp','LineWidth',2)
hold on;
plot(rho, NbIter_6x/Trials,'-m+','LineWidth',2)
plot(rho, NbIter_10x/Trials,'-b*','LineWidth',2)
plot(rho, NbIter_14x/Trials,'-gs','LineWidth',2)
plot(rho, NbIter_x2/Trials,'-rd','LineWidth',2)
xlabel('Sparsity','fontsize',11,'FontWeight','bold');
ylabel('Number of iteration','fontsize',11,'FontWeight','bold');
legend('f(k)=s','f(k)=6k','f(k)=10k','f(k)=14k','f(k)=k^2');
set(gca,'fontsize',11,'FontWeight','bold');
title('Gaussian sparse vectors');
BMPName = 'Number of iteration';
print(gcf,'-dbmp',BMPName)
