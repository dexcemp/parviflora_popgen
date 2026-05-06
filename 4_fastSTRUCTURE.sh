#run fastSTRUCTURE from K: 1 to 10
structure.py -K 1 --input="$prefix"_pruned --output="$prefix"
structure.py -K 2 --input="$prefix"_pruned --output="$prefix"
structure.py -K 3 --input="$prefix"_pruned --output="$prefix"
structure.py -K 4 --input="$prefix"_pruned --output="$prefix"
structure.py -K 5 --input="$prefix"_pruned --output="$prefix"
structure.py -K 6 --input="$prefix"_pruned --output="$prefix"
structure.py -K 7 --input="$prefix"_pruned --output="$prefix"
structure.py -K 8 --input="$prefix"_pruned --output="$prefix"
structure.py -K 9 --input="$prefix"_pruned --output="$prefix"
structure.py -K 10 --input="$prefix"_pruned --output="$prefix"

#finding maximum likelihood K (finding the best K)
python /home/dexcemp/miniconda3/envs/fastStructure/bin/chooseK.py --input="prefix"

#visualization using DISTRUCT
#use the best K results from above
python /home/dexcemp/miniconda3/envs/fastStructure/bin/distruct.py -K 2 --input="prefix" --output="prefix".svg 
