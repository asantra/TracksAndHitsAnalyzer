#! /bin/bash

filename=""

if [[ $1 == "W" ]]
then
   filename="list_root_hics_165gev_w0_3000nm_WIS.txt"
else
   #filename="/nfs/dust/ilc/user/oborysov/hics_list/list_root_hics_165gev_w0_3000nm.txt"
   filename="/nfs/dust/luxe/user/oborysov/hics_list/list_root_hics_165gev_w0_6500nm_jeti40_122020_9550dac4.txt"
fi

# # ### For Tracks tree:
# SECONDS=0
# echo "w0_3000nm JETI40"
# root -l -b -q process_lxtrees_v3.C\(\"${filename}\"\)
# duration=$SECONDS
# echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


#### hits tree
SECONDS=0
echo "w0_5000nm JETI40 Hits"
counter=0
increment=1
rm *HitsProperInfo.txt
while IFS= read -r line
do
  counter=$(( $counter + $increment ))
  echo "$line  $counter"
  echo "$line">oneFileName.txt
  root -l -b -q process_lxtrees_v4.C\(\"oneFileName.txt\",\"$counter\"\)
  break
#   bash eachFileRun.sh oneFileName.txt $counter
done < "$filename"
duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
