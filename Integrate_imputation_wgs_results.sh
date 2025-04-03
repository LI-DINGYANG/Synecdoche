#!/usr/bin/bash

tool=$1
outpath=$2

outdir_Root=${outpath}
outpath=${outpath}/${tool}

# output diretory
if [ ! -d $outpath/ASA ]
then mkdir -p $outpath/ASA
fi

if [ ! -d $outpath/INFO ]
then mkdir -p $outpath/INFO
fi

process_data() {
    local outdir=$1
    local serve=$2
    local i=$3
    local outdir_root=$4
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n' $outdir/chr${i}.dose.vcf.gz | awk -F'\t' '{print $1,$2,$3,$4,$5,$6}' | awk '{if(substr($1, 1, 1) != "c") $1 = "chr"$1; print}' | awk '{if(substr($3, 1, 3) != "chr") $3 = $1":"$2":"$4":"$5; print}' > $outdir/INFO/chr${i}.${serve}.INFO.ASA.txt
    awk '{if ($6 ~ /TYPED(_ONLY)?/) {r2=1} else {split($6,a,";"); for(i=1;i<=length(a);i++) {if(a[i]~/^R2=/) {r2=substr(a[i],4)}}} print $1,$2,$3,$4,$5,r2}' $outdir/INFO/chr${i}.${serve}.INFO.ASA.txt > $outdir/INFO/chr${i}.${serve}.R2.ASA.txt
    rm $outdir/INFO/chr${i}.${serve}.INFO.ASA.txt
    awk 'FNR==NR{a[$1]=$2;next}{if($3 in a){$7=a[$3]} else {$7=-1}}1' $outdir_root/Chinese.WGS1263.frq $outdir/INFO/chr${i}.${serve}.R2.ASA.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.ASA.txt
    bcftools query -f '[%GT ]\n' $outdir/chr${i}.dose.vcf.gz | awk '{for(i=NF-1262;i<=NF;i++) printf("%s ", $i); printf("\n")}' | awk '{for(i=1;i<=NF;i++) gsub(/\|/,"/",$i)} 1' > $outdir/INFO/chr${i}.${serve}.WGSSample.ASA.txt 
    bcftools query -f '[%GP ]\n' $outdir/chr${i}.dose.vcf.gz | awk '{for(i=NF-1262;i<=NF;i++) printf("%s ", $i); printf("\n")}' | awk '{for (i=1; i<=NF; i++) {max=0;split($i, arr, ",");for (j=1; j<=length(arr); j++) {if (arr[j] > max) {max = arr[j];}}$i = max;}print $0;}' > $outdir/INFO/chr${i}.${serve}.IQS.ASA.txt
    paste -d " " $outdir/INFO/chr${i}.${serve}.R2.MAF.ASA.txt $outdir/INFO/chr${i}.${serve}.WGSSample.ASA.txt $outdir/INFO/chr${i}.${serve}.IQS.ASA.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.ASA.txt && rm $outdir/INFO/chr${i}.${serve}.WGSSample.ASA.txt && rm $outdir/INFO/chr${i}.${serve}.IQS.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt
    awk '{for(i=NF-1262-1263;i<=NF-1263;i++) {if($i=="0/0" || $i=="0|0") $i=0; else if($i=="1/1" || $i=="1|1") $i=2; else if($i=="1/0" || $i=="0/1" || $i=="1|0" || $i=="0|1") $i=1; else if($i=="./." || $i==".|.") $i=-1;} print}' $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt
    awk 'NR==FNR{file2_values[$1]=1; next} !($3 in file2_values)' $outdir_root/WGSASA.snplist $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt
    awk 'FNR==NR{a[$3]=$0;next} $3 in a{output = a[$3]; for(i=6; i<=NF; i++) output = output " " $i; print output}' $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt $outdir_root/Chinese.WGS1263.modify.new.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt
    awk '{for(i=8;i<=1270;i++) $i=$i":"$(i+1263+1263); for(i=NF-1262;i<=NF;i++) $i=""; print}' $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt
    awk -v IOSfile="$outdir/INFO/chr${i}.${serve}.IQS.txt" -v labels="0:0 0:1 0:2 0:-1 1:0 1:1 1:2 1:-1 2:0 2:1 2:2 2:-1" '
    BEGIN {
    split(labels, labels_arr, " ")
    for (i = 0; i <= 11; i++) {
        label_probabilities[labels_arr[i]] = 0
    }
    }
    {   
    M =0 
    for (i = 8; i <= NF-1263; i++) {
        label = $i
        prob = $(i + 1263)
        label_probabilities[label] += prob
        if (index(label, "-1") > 0) {
            M += 1;
        }
    }

    N1 = label_probabilities["0:0"] + label_probabilities["1:0"] + label_probabilities["2:0"]
    N2 = label_probabilities["0:1"] + label_probabilities["1:1"] + label_probabilities["2:1"]
    N3 = label_probabilities["0:2"] + label_probabilities["1:2"] + label_probabilities["2:2"]
    S1 = label_probabilities["0:0"] + label_probabilities["0:1"] + label_probabilities["0:2"]
    S2 = label_probabilities["1:0"] + label_probabilities["1:1"] + label_probabilities["1:2"]
    S3 = label_probabilities["2:0"] + label_probabilities["2:1"] + label_probabilities["2:2"]
    N = 1263 - M
    M = 0
    PO = (label_probabilities["0:0"] + label_probabilities["1:1"] + label_probabilities["2:2"]) / N
    PC = (N1 * S1 + N2 * S2 + N3 * S3) / (N*N)
    if (PO == 1 && PC == 1) {
        IQS = 1;
    } else {
        IQS = (PO - PC) / (1 - PC);
    }
    for (i = 0; i <= 11; i++) {
        label_probabilities[labels_arr[i]] = 0
    }
    print IQS,PC >> IOSfile
    for (i = 1; i <= NF - 1263; i++) {
        printf "%s ", $i
    }
    printf "\n"}' $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt

    awk '{for(i=NF-1262; i<=NF; i++){ split($i,a,":"); b[i-NF+1263]=$i" "a[2]" "((a[1]==a[2])? "1" : "0") } output = $1" "$2" "$3" "$4" "$5" "$6" "$7; for(j=1; j<=1263; j++) output = output" "b[j]; print output}' $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt | awk -v x="${serve}" '{$(NF+1)=x; $(NF+2)="ASA"; print}' > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt
    paste -d " " $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt $outdir/INFO/chr${i}.${serve}.IQS.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.ASA.txt && rm $outdir/INFO/chr${i}.${serve}.IQS.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.IQS.txt
    awk '{
    countNotMinusOne = 0
    countCorrect = 0
    for (i = 9; i <= NF-4; i += 3) {
        if ($i!= -1) {
            countNotMinusOne++
            if ($(i + 1) == 1) {
                countCorrect++
            }
        }
    }
    if (countNotMinusOne > 0) {
        ratio = countCorrect / countNotMinusOne
    } else {
        ratio = 0
    }
    print $0, ratio
}' $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.IQS.txt > $outdir/INFO/temp.chr${i}.${serve}.txt && rm $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.IQS.txt && mv $outdir/INFO/temp.chr${i}.${serve}.txt $outdir/INFO/chr${i}.${serve}.R2.MAF.Sample.IQS.txt
}

for chr in {1..22}
do
process_data $outpath $tool $chr $outdir_Root &
done 
wait
#cat $outpath/INFO/*.R2.MAF.Sample.IQS.txt > $outpath/INFO/chr.${tool}.R2.MAF.Sample.IQS.txt


