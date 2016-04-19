% 18/3/2016 
% Trying with modified equations for folding

clc
clear all
close all

PULSAR_PERIOD = 0.03369945724;
% PULSAR_PERIOD = 0.00575731426;
SAMPLING_PERIOD =0.0001;
SAMPLESIZE =1536;
NSAMPLES =10000;
NBINS  =input('Enter the number of bins. \n');
NFILES =input('Enter the number of files. \n');
% i = bin number
% j = sample number
% for NBINS = NOBINS:50+NOBINS
B = zeros(SAMPLESIZE,10000);
C = zeros(SAMPLESIZE,NBINS);
D = zeros(SAMPLESIZE,NBINS);
E = zeros(SAMPLESIZE,NBINS); 
counter =zeros(NBINS,1);

for p=7068:7067+NFILES
    
    filename=strcat('20130924B/20130924B.137996',num2str(p),'.beam');
    fileID1 =fopen(filename);
    A = fread(fileID1,15360000,'uint8=>uint8');
    fclose(fileID1);
    
    k=0;
    l=0;
    t=0;
 for j=1:10000
 B(:,j) =double(A(((j-1)*SAMPLESIZE)+1:(j*SAMPLESIZE),1));
 end
 
 for j=1:1536
    a=sum(B(j,:))/10000;
    for n=1:10000
        C(j,n)=B(j,n)-a;
    end 
    a=0;
end
    t =ceil(PULSAR_PERIOD*NSAMPLES);
    for j=1:10000
        j;
        m =ceil(j/t);
        k =(j-1)*SAMPLING_PERIOD/PULSAR_PERIOD-(m-1);
        if k>1 
          k=0;
        end
        i = floor (k*NBINS)+1;
        D(:,i) = D(:,i)+ C(:,j);
        counter(i,:) = counter(i,:)+1;
    end
    
for j=1:NBINS
    D(:,j)= D(:,j)/counter(j,:);
end
   
    
for j=1:1536
    a=sum(D(j,:))/NBINS;
    for n=1:NBINS
        E(j,n)=D(j,n)-a;
    end 
    a=0;
end
end
% counter(:,:)
% sum(counter(:,:))

nu=1:NBINS;
t=1:1536;
figure
pcolor(nu,t,E);
shading interp
colorbar

% end