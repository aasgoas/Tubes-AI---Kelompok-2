## Mutasi menggunakan swapping Mutasi

## Input:
## Kromosom :kromosom, matriks berukuran 1x JumGen
## JumGen: jumlah gen
## Pmutasi: probabilitas mutasi

## Ouput:
## MutKrom : kromosom hasil mutasi, matriks berukuran 1 x JumGen

function MutKrom = TSPMutasi(Kromosom,JumGen,Pmutasi)
MutKrom =Kromosom;
for ii=1:JumGen,
if rand <Pmutasi,
45
 TM2=1+fix(rand*JumGen);
while TM2==ii,
 TM2=1+fix(rand*JumGen);
end
 temp = MutKrom(ii);
 MutKrom(ii)=MutKrom(TM2);
 MutKrom(TM2)=temp;
end
end
