clear all
 close all


load('SS_N1.mat');

%% load SS signals, tsec time, NN number of signals

%z-normalize:

MeanSS= mean(SS,2);
StdSS=std(SS,1,2);
for i=1:1:NN
zSS(i,:)=(SS(i,:)-MeanSS(i))/StdSS(i);
end


% Amplitude istogram:

nbin=101;
xmin=-5;
xmax=5;
passo=(xmax-xmin)/nbin;

b=(xmin+passo/2):passo:(xmax-passo/2);

Ndat=size(zSS,2);

hmedia=zeros(size(b));

figure
clear zSSk0;
med=0;
for k=1:1:NN
    clear HH
    clear h
    med=med+1;
     zSSk0=zSS(k,:);

HH=zeros(size(b));
    for i=1:1:Ndat
    h=1+floor(nbin*(zSSk0(i)-xmin)/(xmax-xmin));
    if h>=1 && h<=nbin
        HH(h)=HH(h)+1;
    end
    end

    for h=1:1:nbin
    HH(h)=HH(h)*nbin/(Ndat*(xmax-xmin));  %%% normalizzazione tale che integrale fa uno
    hmedia(h)=hmedia(h)+HH(h);
    end

    hold on
  
plot(b,HH,'.b')
    
    end
   
   legend('show') 


 
hmedia=hmedia/med;
plot(b,hmedia,'b-','LineWidth',2,'Displayname','average of data')

hold on

% comparison with gaussian 

mu=0;
sigma=1/(hmedia(int8(nbin/2))*sqrt(2*pi));
%calcolo sigma per avere stessa altezza dei dati medi 
hold on
y=-5:0.1:5;

cc=1;
f=cc*exp(-(y-mu).^2/(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f,'r-','LineWidth',2,'Displayname','GAUSSIAN ')


%% FROM amplitude flucutation picture choose threshold

SOGLIA_AMP=2;  %threshold
fprintf("threshold  %f \n",SOGLIA_AMP)

for i=1:NN

    clear select_zSS;
select_zSS=abs(zSS(i,:))-SOGLIA_AMP;
% sottraggo SOGLIA_AMP e metto a zero quella parte che viene negativa perche non superava SOGLIA_AMP !!!
select_zSS(select_zSS <=0 )=0;

% prendo il valore assoluto dei segnali che superano la soglia
selezionato_zSS(i,:)=abs(select_zSS(:));
end

somma=sum(abs(selezionato_zSS),1);

clear valanga
clear valanga2
clear vA;
clear taglia;
soglia_somma=0;

k=0;
startvalanga=0;
for j=2:1:size(tsec,2)
      if (somma(j)>soglia_somma & somma(j-1)<= soglia_somma  )      
        startvalanga=tsec(j);
        startj=j;
         k=k+1 ;
      elseif (somma(j)<= soglia_somma & somma(j-1)> soglia_somma)
         if startvalanga>0 
         valanga(k)=tsec(j)-startvalanga;
         clear pezzosom;
         pezzosom=somma(startj:j);
          taglia(k)=sum(pezzosom);
         end
       end
end
 
% avalanche's size is written in: taglia
%avalanche's duration is written in: valanga

flagValanga=zeros(size(tsec));
flagValanga(somma>soglia_somma)=1;




%%% HISTO SIZE

dur1=taglia';  %vector with all avalanches' sizes.


fprintf("number of avalanches ? %i \n",size(dur1,1))


% histo lineare

xmin=0;
xmax=max(dur1);
nbin=35000;
passo=(xmax-xmin)/nbin;

histo1=zeros(nbin,1);
xx=zeros(nbin,1);

Ndat=size(dur1,1);
for i=1:1:Ndat;
    h=ceil(nbin*(dur1(i)-xmin)/(xmax-xmin));
    if h>=1 && h<=nbin;
        histo1(h)=histo1(h)+1;
    end
end

for i=1:1:nbin;
    xx(i)=xmin+(i-0.5)*(xmax-xmin)/nbin;
    histo1(i)=histo1(i)*nbin/(Ndat*(xmax-xmin));
end

histo=histo1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear b;
xmin=0.1;
xmax=max(dur1);
nbin=30;
%nbin=17
b=1:1:nbin;
norml=zeros(size(b));

for i=1:nbin;
    b(i)=xmin*exp((i-0.5)*log(xmax/xmin)/(nbin));
    norml(i)=xmin*(exp(i*log(xmax/xmin)/(nbin))-exp((i-1)*log(xmax/xmin)/(nbin)));
end

HH=zeros(size(b));
Ndat=size(dur1,1);
for i=1:1:Ndat;
    h=ceil(nbin*log(dur1(i)/xmin)/log(xmax/xmin));
    if h>= 1 && h<= nbin;
        HH(h)=HH(h)+1;
    end
end

for i=1:1:nbin;
    HHS(i)=(HH(i))/(norml(i)*Ndat);   % normalizzazione
end

bS=b;

xxS=xx;
histoS=histo;

figure

loglog(xx,histo,'*','color','#4DBEEE','MarkerSize',15,'LineWidth',0.1)
hold on
loglog(bS,HHS,'ro','MarkerSize',13,'LineWidth',3)
hold on

xlabel('S ')
ylabel('P(S)  ')
hold on
legend( 'histo linear binning ',' histo log-bin' )

axis([0.8 10000 0.000000008 1])

%%%%%
%%%%% HISTO  dutations T

clear h*
clear p*
clear b*
clear d*
clear x*
clear n*

dur1=valanga';

xmin=0;
xmax=max(dur1);
nbin=3000 % number of LINEAR BINS
passo=(xmax-xmin)/nbin;

histo1=zeros(nbin,1);
xx=zeros(nbin,1);

Ndat=size(dur1,1);
for i=1:1:Ndat;
    h=ceil(nbin*(dur1(i)-xmin)/(xmax-xmin));
    if h>=1 && h<=nbin;
        histo1(h)=histo1(h)+1;
    end
end

for i=1:1:nbin;
    xx(i)=xmin+(i-0.5)*(xmax-xmin)/nbin;
    histo1(i)=histo1(i)*nbin/(Ndat*(xmax-xmin));
end
histo=histo1;

clear b;

xmin=(1/128);
xmax=max(dur1);
nbin=30; % number of Logaritmic BINS

b=1:1:nbin;
norml=zeros(size(b));

for i=1:nbin;
    b(i)=xmin*exp((i-0.5)*log(xmax/xmin)/(nbin));
    norml(i)=xmin*(exp(i*log(xmax/xmin)/(nbin))-exp((i-1)*log(xmax/xmin)/(nbin)));
end

HH=zeros(size(b));
Ndat=size(dur1,1);
for i=1:1:Ndat;
    h=ceil(nbin*log(dur1(i)/xmin)/log(xmax/xmin));
    if h>= 1 && h<= nbin;
        HH(h)=HH(h)+1;
    end
end

for i=1:1:nbin;
    HHS(i)=(HH(i))/(norml(i)*Ndat);   % normalizzazione
end

bS=b;

figure  % figure P(T)

loglog(xx,histo,'*','color','#4DBEEE','MarkerSize',15,'LineWidth',1)
hold on
loglog(bS,HHS,'ro','MarkerSize',13,'LineWidth',3)
hold on

xlabel('T ')
ylabel('P(T)  ')
legend( 'histo linear binning','logaritmic binning' )
axis([0.024 20.1 0.00001 100])

