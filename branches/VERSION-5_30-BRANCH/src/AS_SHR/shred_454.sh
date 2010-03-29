#!/bin/sh

# Take in an ACE file generated by the 454 assembly software.
# Shred the contigs into pseudo reads, 600bp long with 300bp overlaps.
# Write an FRG file for intput to Celera Assembler.

ACE=$1
echo "Input ACE file was..."
echo $ACE
ls -l $ACE

if test $2
then
QUAL=$2
else
QUAL=3
fi
echo "Quality score is set to be"
echo $QUAL


BINDIR=/usr/local/packages/CA/bin

echo ">>Generate_NonShallow_Contigs..."
$BINDIR/Generate_NonShallow_Contigs.pl -a $ACE -f pyro.contigs

echo ">>Shred_Contigs..."
$BINDIR/Shred_Contigs.pl -f pyro.contigs > pyro.shreds

echo ">>FASTA_to_frg_file..."
$BINDIR/FASTA_to_frg_file.pl -f pyro.shreds  -q $QUAL> pyro.frg

echo ">>Done"