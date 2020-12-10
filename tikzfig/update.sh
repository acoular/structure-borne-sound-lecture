#!/bin/bash
for pdf in *.pdf
do
svg="${pdf%.*}".svg
if ! [ -e $svg  ] || [ $(($(date -r $svg +%s)-$(date -r $pdf +%s) )) -le 0 ]; then 
echo $pdf
pdf2svg $pdf $svg
fi
done
