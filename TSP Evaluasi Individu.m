## Mengevaluasi sehingga mendapatkan nilai fitness nya
## Input:
## Kromosom: kromosom dengan matrik ukuran 1XJumGen
## XYKota : posisi-posisi koordinat kota, matrik berukuran UkPop x 2
## JumGen: jumlah gen
## Output:
## fitnnes : nilai fitness

function Fitness = TSPEvaluasiIndividu(Kromosom,JumGen,XYkota)
TB = 0;
for ii=1:JumGen-1,
 TB=TB+norm(XYkota(Kromosom(ii),:)-XYkota(Kromosom(ii+1),:));
end
% Jalur untuk kembali ke kota asal.
TB =TB +norm(XYkota(Kromosom(ii),:)-XYkota(Kromosom(1),:));
Fitness =1/TB;
