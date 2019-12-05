#!/bin/bash
#SBATCH --job-name=ErrorRateCalc
#SBATCH --time=96:00:00
#SBATCH -p normal

module load parallel
module load bio-misc
module load R/gcc7/3.6.1

ErrorRateCalc () {
fqBase=$1
for barcode in CGATGCTCTGCA AAGCCGGTTGCA; do
        for read in R1 R2; do
        agrep -2 -D100 -I100 "^${barcode}" ${fqBase}.${read}.fq | head -n -1 | sed 's/^/2del,--/g' > ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^.${barcode}" ${fqBase}.${read}.fq | head -n -1 | sed 's/^/1del,-/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^..${barcode}" ${fqBase}.${read}.fq | head -n -1 | sed 's/^/0Ind,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^...${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 2- | sed 's/^/1Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^....${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 3- | sed 's/^/2Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^.....${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 4- | sed 's/^/3Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^......${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 5- | sed 's/^/4Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^.......${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 6- | sed 's/^/5Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^........${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 7- | sed 's/^/6Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^.........${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 8- | sed 's/^/7Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq
        agrep -2 -D100 -I100 "^..........${barcode}" ${fqBase}.${read}.fq | head -n -1 | cut -c 9- | sed 's/^/8Ins,/g' >> ${fqBase}.${read}.${barcode}.NoDels.seq

        cut -c 1-21 ${fqBase}.${read}.${barcode}.NoDels.seq | sed -e 's/\([GATCN-]\)/\1,/g' | sed -e 's/,$//g' -e "s/^/$fqBase,$barcode,$read,/g" > ${fqBase}_${barcode}_${read}.csv #save to file

        #Output Summary Stats
        # echo ${fqBase},${read},${barcode},$P1del,$P2del,$PropP1Del,$PropP2Del,$PropNoDels,$PropOneIns,$PropTwoIns >> ${fqBase}_${barcode}_${read}_Summary.csv #add more
        done
done

}
export -f ErrorRateCalc
#run function in parallel
ls *R1.fq | sed 's/\.R1\.fq//g' | parallel -j 7 --no-notice ErrorRateCalc {}

Albatross_Barcodes=$(ls *Albatross*[1-2].csv)
Contemp_Barcodes=$(ls *Contemp*[1-2].csv)

Header=fqBase,Barcode,Read,Indels,BP1,BP2,BP3,BP4,BP5,BP6,BP7,BP8,BP9,BP10,BP11,BP12,BP13,BP14,BP15,BP16

cat <(echo $Header) $Albatross_Barcodes > AlbatrossBarcodes.csv
cat <(echo $Header) $Contemp_Barcodes > ContempBarcodes.csv

R CMD BATCH PIRE_Aduo_Stats_forR_12_2_19.R
