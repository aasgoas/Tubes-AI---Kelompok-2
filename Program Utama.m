clear
clc
XYkota = [1 3; 1 7; 3 9; 5 3; 7 1; 9 5];
JumGen =length(XYkota(:,1)); 
%%Jumlah Gen
UkPop =100;
%%Jumlah Kromosom dalam Populasi
Psilang=0.8;
%%Probabilitas Pindah Silang
Pmutasi=0.005;
%%Probabilitas Mutasi
MaxG =100;
%%Jumlah Generasi
PanjJalHarp = 40;
Fthreshold = 1/PanjJalHarp;
Bgraf=Fthreshold;
%%Inisialisasi Grafis
hfig=figure;
hold on
set(hfig,'position',[50,50,600,400]);
set(hfig,'DoubleBuffer','on');
axis([1 MaxG 0 Bgraf]);
hbestplot1 = plot(1:MaxG,zeros(1,MaxG));
hbestplot2 = plot(1:MaxG,zeros(1,MaxG));
htext1 =text(0.6*MaxG,0.25*Bgraf, sprintf('Fitness terbaik: %7.6f', 0.0));
htext2 =text(0.6*MaxG,0.20*Bgraf, sprintf('Fitness rata-rata: %7.6f', 0.0));
htext3 =text(0.6*MaxG,0.15*Bgraf, sprintf('Panjang Jalur Terbaik: %7.3f', 0.0));
htext4 =text(0.6*MaxG,0.10*Bgraf, sprintf('Ukuran Populasi: %3.0f', 0.0));
htext5 =text(0.6*MaxG,0.05*Bgraf, sprintf('Probabilitas Mutasi: %4.3f', 0.0));
xlabel('Generasi');
ylabel('Fitness');
hold off
drawnow;

%%Inisialisasi Pupulasi
function Populasi = TSPInisialisasiPopulasi(UkPop,JumGen)
for ii=1:UkPop,
        [Xval,Ind]= sort(rand(1, JumGen));
        Populasi(ii,:) =Ind;
end
for generasi=1:MaxG,
        MaxF= TSPEvaluasiIndividu(Populasi(1,:),JumGen,XYkota);
        MinF = MaxF;
        IndeksIndividuTerbaik=1;
for ii=1:UkPop,
    Fitness(ii) = TSPEvaluasiIndividu(Populasi(ii,:),JumGen,XYkota);
if (Fitness(ii)> MaxF),
            MaxF = Fitness(ii);
            IndeksIndividuTerbaik=ii;
            JalurTerbaik = Populasi(ii,:);
end
if (Fitness(ii)<= MinF)
            MinF = Fitness(ii);
end
end
end
    FitnessRata = mean(Fitness);
    plotvector1 = get(hbestplot1,'YData');
    plotvector1(generasi)= MaxF;
    set(hbestplot1,'YData',plotvector1);
    plotvector2 = get(hbestplot2,'YData');
    plotvector2(generasi)= FitnessRata;
    set(hbestplot2,'YData',plotvector2);
    set(htext1,'String',sprintf('Fitness terbaik: %7.6f',MaxF));
    set(htext2,'String',sprintf('Fitness rata-rata: %7.6f',FitnessRata));
    set(htext3,'String',sprintf('Panjang Jalur Terbaik: %7.3f',1/MaxF));
    set(htext4,'String',sprintf('Ukuran Populasi: %3.0f',UkPop));
    set(htext5,'String',sprintf('Probabilitas Mutasi: %4.3f',Pmutasi));
    drawnow;
if MaxF >Fthreshold,
RETURN;
end
    TemPopulasi = Populasi;
%%Elitisme
%%-Buat satu kopi kromosom terbaik jika ukuran populasi ganjil
%%-Buat 2 kopi kromosom terbaik jika ukuran populasi genap
if mod(UkPop,2)==0,
        IterasiMulai =3;
        TemPopulasi(1,:) = Populasi(IndeksIndividuTerbaik,:);
        TemPopulasi(2,:) = Populasi(IndeksIndividuTerbaik,:);
else
        IterasiMulai =2;
        TemPopulasi(1,:) = Populasi(IndeksIndividuTerbaik,:);
end
 LinearFitness =LinearFitnessRanking(UkPop,Fitness,MaxF,MinF);

%%Roulette-wheel dan Pindah Silang
for jj=IterasiMulai:2:UkPop,
        IP1 =RouletteWheel(UkPop,LinearFitness);
        IP2 =RouletteWheel(UkPop,LinearFitness);
if (rand<Psilang),
        Anak = TSPPindahSilang(Populasi(IP1,:),Populasi(IP2,:),JumGen);
        TemPopulasi(jj,:) = Anak(1,:);
        TemPopulasi(jj+1,:) = Anak(2,:);
else
        TemPopulasi(jj,:) = Populasi(IP1,:);
        TemPopulasi(jj+1,:) = Populasi(IP2,:);
end
end
%%Mutasi pada keseluruhan kromosom
for kk=IterasiMulai:UkPop,
    TemPopulasi(kk,:) = TSPMutasi(TemPopulasi(kk,:), JumGen, Pmutasi);
end
    Populasi =TemPopulasi;
    save JalurTerbaik.matJalurTerbaik
end
%%Tanpa tanda ';' berarti menampilkan nilai dari variabel 
