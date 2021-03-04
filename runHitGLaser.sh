#! /bin/bash


SECONDS=0

filename=""

if [[ $1 == "W" ]]
then
    filename="gPlusLaserWIS.txt"
else
#    filename="gPlusLaser.txt"
#     filename="/nfs/dust/luxe/user/oborysov/hics_list/list_root_bppp_background_fast_0508546b_0_5.txt"
    filename="/nfs/dust/luxe/group/MCProduction/Background/gammalaser/09102020_lxb18e/Merged/Files"
fi

# # For Tracks tree:
# echo "photon+laser background samples"
# root -l -b -q process_track_tree_draw_v5.C\(\"${filename}\"\)
# duration=$SECONDS
# echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
# 
# exit 1
#### hits tree
SECONDS=0
echo "photon+laser background samples"
# root -l -b -q process_lxtrees_v2.C\(\"${filename}\"\)
root -l -b -q process_lxtrees_background.C\(\"${filename}\"\)
duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
