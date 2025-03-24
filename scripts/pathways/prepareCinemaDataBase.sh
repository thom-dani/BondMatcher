#!/bin/sh

DataBaseFile=$1

if [ -z "$DataBaseFile" ]; then
  echo "Usage:"
  echo "  ./prepareCinemaDataBase.sh <database file>"
  exit 0
fi

echo "Decompressing archive $DataBaseFile..."
tar xvzf $DataBaseFile 

DataBase=${DataBaseFile/.tar.gz/.cdb}

echo "Creating Cinema database $DataBase..."

mkdir $DataBase
mv data $DataBase/
mv helper.txt $DataBase/data.csv
