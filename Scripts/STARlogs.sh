#!/bin/bash

mkdir Logs;

for DIR in `ls STAR`;
do 
	echo -ne "$DIR\t" >> Logs/Input.tsv;
	cat STAR/$DIR/Log.final.out | grep "Number of input reads" | cut -f2 >> Logs/Input.tsv;
done

for DIR in `ls STAR`;
do 
	echo -ne "$DIR\t" >> Logs/Unique.tsv;
	cat STAR/$DIR/Log.final.out | grep "Uniquely mapped reads number" | cut -f2 >> Logs/Unique.tsv;
done

for DIR in `ls STAR`;
do 
	echo -ne "$DIR\t" >> Logs/Multiple.tsv;
	cat STAR/$DIR/Log.final.out | grep "Number of reads mapped to multiple loci" | cut -f2 >> Logs/Multiple.tsv;
done

for DIR in `ls STAR`;
do 
	echo -ne "$DIR\t" >> Logs/TooMany.tsv;
	cat STAR/$DIR/Log.final.out | grep "Number of reads mapped to too many loci" | cut -f2 >> Logs/TooMany.tsv;
done

join Logs/Input.tsv Logs/Unique.tsv | join - Logs/Multiple.tsv | join - Logs/TooMany.tsv > Logs/Logs.tsv


