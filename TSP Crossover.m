## Pindah Silang mengggunakan Order Crossover

## Input:
## Bapak :kromosom, matriks berukuran 1x JumGen
## Ibu :kromosom, matriks berukuran 1x JumGen
## JumGen: jumlah gen

## Output:
## Anak : kromosom hasil pindah silang, matriks berukuran 1 x 

function Anak = TSPPindahSilang(Bapak,Ibu,JumGen)
cp1 = 1 + fix(rand*(JumGen-1));
cp2 = 1 + fix(rand*(JumGen-1));
while cp2==cp1,
 cp2=1+fix(rand*(JumGen-1));
end
if cp1<cp2,
 cps =cp1;
 cpd =cp2;
44
else
 cps =cp2;
 cpd =cp1;
end
Anak(1,cps+1:cpd)=Ibu(cps+1:cpd);
Anak(2,cps+1:cpd)=Bapak(cps+1:cpd);
SisaGenbapak =[];
SisaGenibu =[];
for ii=1:JumGen,
if ~ismember(Bapak(ii), Anak(1,:));
 SisaGenbapak = [SisaGenbapak Bapak(ii)];
end
if ~ismember(Ibu(ii), Anak(1,:));
 SisaGenibu = [SisaGenibu Ibu(ii)];
end
end
Anak(1,cpd+1:JumGen) = SisaGenbapak(1:JumGen-cpd);
Anak(1,1:cps) = SisaGenbapak(1+JumGen-cpd:length(SisaGenbapak));
Anak(2,cpd+1:JumGen) = SisaGenibu(1:JumGen-cpd);
Anak(2,1:cps) = SisaGenibu(1+JumGen-cpd:length(SisaGenibu));
