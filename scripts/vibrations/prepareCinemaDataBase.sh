#!/bin/sh

DataBaseFile1=$1
DataBaseFile2=$2
DataBaseFile3=$3
DataBaseFile4=$4

if [ -z "$DataBaseFile4" ]; then
  echo "Usage:"
  echo "  ./prepareCinemaDataBase.sh <file1> <file2> <file3> <file4>"
  exit 0
fi

DataBaseFile=${DataBaseFile1/.partaa/}

echo "Merging files into $DataBaseFile..."
cat $DataBaseFile1 $DataBaseFile2 $DataBaseFile3 $DataBaseFile4 >  $DataBaseFile

echo "Decompressing archive $DataBaseFile..."
tar xvzf $DataBaseFile 

DataBase=${DataBaseFile/.tar.gz/.cdb}

echo "Creating Cinema database $DataBase..."

mkdir $DataBase
mv data $DataBase/
mv helper.txt $DataBase/data.csv
