#!/bin/sh

if [ "$#" -lt "2" ]; then 
  echo "Usage: ./separatrixExtraction.sh <input db> <resolution> [<save preprocessed data, -s>]"
  exit 0
fi

input=$1

resolution=$2

savePreprocessedData="0"
if [ "$3" == "-s" ]; then
  savePreprocessedData="1"
fi

coreNumber=`nproc`

#coreNumber=`echo "$coreNumber / 2" | bc`

output=${input}_resolution$resolution.cdb/
if [ "${input: -1}" == "/" ]; then
  # finishes with a slash
  output=${input::-1}_resolution$resolution.cdb/
fi

echo "--------------------------------------------------------------------------------"
echo "Runnning separatrix extraction:"
echo "  Input: $input"
echo "  Output: $output"
echo "  Resolution: $resolution"
echo "  Core number: $coreNumber"
echo "  Save preprocessed data: $savePreprocessedData"
echo "--------------------------------------------------------------------------------"

if [ ! -e $output ]; then

  echo -e "\n\n\nCopying directory structure for $input into $output..."
  find $input -type d | xargs -I{} mkdir -p "_tmp/{}"
  mv _tmp/$input $output
  rm -R _tmp

  echo -e "\n\n\nCopying database index to $output..."
  cp $input/data.csv $output/data.csv
  sed -i 's/start_data.vti/start_data.vtu/g'  $output/data.csv

fi

initDir=`pwd`

cd $input

for i in data/*/*/*/*/*/*/*/*/*/*; do
  if [ ! -e "$initDir/$output/${i}/start_data.vtu" ]; then
    echo -e "\n\n\nProcessing entry ${i}..."

    inputFile=$initDir/$input/$i/start_data.vti

    cd $initDir/$output/$i
    ln -sf $initDir/separatrixExtraction.py .
    ln -sf $initDir/separatrixPreprocess.py .
    ln -sf $inputFile

    if [ ! -e "$initDir/$input/$i/input${resolution}.pvti" ]; then
      echo -e "\tRunning pre-processing..."
      mpirun --oversubscribe -n $coreNumber pvbatch separatrixPreprocess.py $resolution 
      mv input${resolution} $initDir/$input/$i/input${resolution}
      mv input${resolution}.pvti $initDir/$input/$i/input${resolution}.pvti
      ln -sf $initDir/$input/$i/input${resolution} .
      ln -sf $initDir/$input/$i/input${resolution}.pvti
    fi

    echo -e "\tExtracting separatrices..."
    ./separatrixExtraction.py input${resolution}.pvti  &> separatrixExtraction.log

    if [ "$savePreprocessedData" == "0" ]; then
      rm -R $initDir/$input/$i/input${resolution}
      rm $initDir/$input/$i/input${resolution}.pvti
      rm input${resolution}
      rm input${resolution}.pvti
    fi

    cd $initDir/$input

    # debug
    #cd $initDir
    #exit 0
  fi
done

cd $initDir
