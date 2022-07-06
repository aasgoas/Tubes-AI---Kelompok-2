%%Membangkitkan ukran populasi kromosom masing-masing kromosom 
berisi
%%bilangan 1 sampai JumGen(jumlah kota)
%%
%%Input:
%%UkPop : Ukuran populasi atau jumlah kromosom dalm populasi
%%JumGen: jumlah gen(jumlah kota)
%%
%%Output:
%%Pupulasi : kumpulan kromosom matriks berukuran UkPop x JumGen
%%
function Populasi = TSPInisialisasiPopulasi(UkPop,JumGen)
for ii=1:UkPop,
 [Xval,Ind]= sort(rand(1, JumGen));
 Populasi(ii,:) =Ind;
end
