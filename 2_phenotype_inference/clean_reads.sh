#! /bin/sh -e

clean () (
        gunzip -c $* |
        awk -F : 'BEGIN {OFS="\t"}
                NR%4==1 {idx=$10}
                NR%4==2 {print idx, $0}
        ' |
        tr '+' '\t'
)

mkfifo ltseq-paste-$$.1
mkfifo ltseq-paste-$$.2

clean $1 > ltseq-paste-$$.1 &
clean $2 | awk '{print $3}'  > ltseq-paste-$$.2 &

#clean $1/$2*R1.fastq* > ltseq-paste-$$.1 &
#clean $1/$2*R2.fastq* | awk '{print $3}'  > ltseq-paste-$$.2 &

paste ltseq-paste-$$.1 ltseq-paste-$$.2 |
grep -v -P 'G{10}(      |$)'
wait
rm -rf ltseq-paste-$$.1 ltseq-paste-$$.2
