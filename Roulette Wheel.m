## Memilih orang tua menggunkan Linear Fitness yaitu nilai fitness hasil
## pen-skala-an. Pemilihan dilakukan secara proporsional sesuai dengan nilai
## fitness-nya

## Input:
## UkPop : ukuran populasi atau jumlah kromosom dalam populasi
## LinearFitness : nilai fitness yang sudah di skala-kan

## Output :
## Pindex : indeks dari kromosom yang terpilih

function Pindex = RouletteWheel(UkPop,LinearFitness);
JumFitness = sum(LinearFitness);
KumulatifFitness =0;
RN = rand;
ii = 1;
while ii<=UkPop,
 KumulatifFitness =KumulatifFitness + LinearFitness(ii);
if (KumulatifFitness/JumFitness)> RN,
 Pindex =ii;
break;
end
 ii =ii+1;
end
