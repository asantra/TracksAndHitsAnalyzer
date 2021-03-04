#! /bin/bash



# For Tracks tree:
echo "e+laser electron only samples"

filename=""

if [[ $1 == "W" ]]
then
   #filename="EBeamOnlyWIS.txt"
#    filename1="EBeamOnlyWIS_DividedByBX1.txt"
#    filename2="EBeamOnlyWIS_DividedByBX2.txt"
#    filename3="EBeamOnlyWIS_DividedByBX3.txt"
#    filename4="EBeamOnlyWIS_DividedByBX4.txt"
    filename="QGSPList_EBeamOnly.txt"
else
#    #filename="EBeamOnly.txt"
#    filename="/nfs/dust/luxe/user/oborysov/hics_list/list_root_hics_background_fast_c99bba6d_0_18_luxe.txt"
#     filename="/nfs/dust/luxe/group/MCProduction/Background/elaser/29102020_lx86a1/Merged/Files/"
    filename="/nfs/dust/luxe/user/oborysov/hics_list/list_root_hics_background_qgsp_9a61db54_0_13.txt"
fi

# SECONDS=0
# root -l -b -q process_track_tree_draw_v6.C\(\"${filename}\"\)
# duration=$SECONDS
# echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

SECONDS=0
if [[ $1 == "W" ]]
then
#     root -l -b -q process_lxtrees_v2.C\(\"${filename1}\"\)
#     root -l -b -q process_lxtrees_v2.C\(\"${filename2}\"\)
#     root -l -b -q process_lxtrees_v2.C\(\"${filename3}\"\)
#     root -l -b -q process_lxtrees_v2.C\(\"${filename4}\"\)
    root -l -b -q process_lxtrees_background.C\(\"${filename}\"\)
else
#     root -l -b -q process_lxtrees_v2.C\(\"${filename}\"\)
    root -l -b -q process_lxtrees_background.C\(\"${filename}\"\)
fi

duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
