## Menskalakan nilai fitness ke dalam ranking sehingga diperoleh nilai-nilai
## fitness baru yang berada dalam rentang (MaxF,MinF)

## Input:
## UkPop : ukuran populasi atau jumlah kromosom dalam populasi
## Fitness: nilai fitness, matriks ukuran 1x UkPop
## MaxF : Nilai Fitness Maximum 
## MinF : Nilai Fitness minimum 

## Output :
## LFR : Linear Fitness Ranking

function LFR = LinearFitnessRanking(UkPop,Fitness,MaxF,MinF);
[SF,IndF]=sort(Fitness);
for rr=1; UkPop,
 LFR (IndF(UkPop-rr+1)) = MaxF-(MaxF-MinF)*((rr-1)/(UkPop1));
