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

# dir structure:
# - data
#   - W6_drop1_struc_<pathwayParameter>
#     - nonrel/PBE0/TZ2P/gridspacing0.05/rho/

for i in data/*; do

  outputDir="$initDir/$output/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/"

  if [ ! -e "$outputDir/start_data.vtu" ]; then
    echo -e "\n\n\nProcessing entry ${i}..."
    inputFile="$initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/start_data.vti"
    cd $outputDir
    ln -sf $initDir/separatrixExtraction.py .
    ln -sf $initDir/separatrixPreprocess.py .
    ln -sf $inputFile

    inputPreprocessedFile="$initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/input${resolution}.pvti"

    if [ ! -e "$inputPreprocessedFile=" ]; then
      echo -e "\tRunning pre-processing..."
      mpirun --oversubscribe -n $coreNumber pvbatch separatrixPreprocess.py $resolution
      mv input${resolution} $initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/input${resolution}
      mv input${resolution}.pvti $initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/input${resolution}.pvti
      ln -sf $initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/input${resolution} .
      ln -sf $initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/input${resolution}.pvti
    fi

    echo -e "\tExtracting separatrices..."
    ./separatrixExtraction.py input${resolution}.pvti  &> separatrixExtraction.log

    if [ "$savePreprocessedData" == "0" ]; then
      rm -R $initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/input${resolution}
      rm $initDir/$input/$i/nonrel/PBE0/TZ2P/gridspacing0.05/rho/input${resolution}.pvti
      rm input${resolution}
      rm input${resolution}.pvti
    fi

    cd $initDir/$input
  fi

done

cd $initDir
