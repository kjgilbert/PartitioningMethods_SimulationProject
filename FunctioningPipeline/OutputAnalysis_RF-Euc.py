# input / output
import sys

# euclidean distance calculation
import dendropy
from dendropy import Tree as tree
from dendropy.calculate import treecompare


# Jeremy's RF calculation script
from ete3 import Tree


def collapse_branches(TREE_FILE,SUPPORT_THRESHOLD):
    t = Tree(TREE_FILE)
    for node in t.get_descendants():
        if not node.is_leaf() and (node.support <= SUPPORT_THRESHOLD):
            node.delete()
    return t


def count_internal(tree):
    tree.unroot()
    edges=-1
    for edge in tree.traverse():
        if not edge.is_leaf():
            edges+=1
    return edges

def rf_distance(tree1,tree2,option=False):
    if option=='collapse':
        t1 = collapse_branches(tree1,0.75)
        t2 = collapse_branches(tree2,0.75)
        option = 'reduced'
    else:
        t1 = tree1 #Tree(tree1) because I import the trees from files outside this function
        t2 = tree2 #Tree(tree2)

    t1.unroot()
    t2.unroot()

    rf = t1.robinson_foulds(t2,unrooted_trees=True)

    rf_dist = rf[0]
    max_rf = rf[1]
    num_leaves=len(rf[2])

    max_resolved_score = (2*num_leaves)-6

    internal_branches_t1 = count_internal(t1)
    internal_branches_t2 = count_internal(t2)

    total_internal = internal_branches_t1 + internal_branches_t2

    num_missing_splits_t1 = num_leaves - 3 - internal_branches_t1
    num_missing_splits_t2 = num_leaves - 3 - internal_branches_t2



    rf_dist_upper = rf_dist + num_missing_splits_t1 + num_missing_splits_t2


    #If normalise by (num intenral branches)
    if option=='reduced':
        normalised_rf = rf_dist/total_internal

    elif option=='upper':
    #If add score to create upper bound due to being polytomy
        normalised_rf = rf_dist_upper/max_resolved_score
    else:
        normalised_rf = rf_dist/max_resolved_score

    return normalised_rf

#________________________________________________________________________________________#


## REFORMAT THE IDs IN THE REAL TREE SO FOLLOWING COMPARISONS WILL WORK
## HAPPENS IN BASH RunPipeline_Script.sh


## EXTRACT THE INFERRED TREES FROM THE OUTPUT FILE
## HAPPENS IN BASH RunPipeline_Script.sh


#________________________________________________________________________________________#

# analyses on the outputs, using the above functions:

inputTree1 = sys.argv[1]
inputTree2 = sys.argv[2]

realTree = Tree(inputTree1)
inferredTree = Tree(inputTree2)

# get RF distance
rf_dist = (rf_distance(realTree, inferredTree))



# get Euclidean distance

tns = dendropy.TaxonNamespace()

tree1 = tree.get_from_path(inputTree1, "newick", taxon_namespace=tns)
tree2 = tree.get_from_path(inputTree2, "newick", taxon_namespace=tns)

tree1.encode_bipartitions()
tree2.encode_bipartitions()

eucl_dist = (treecompare.euclidean_distance(tree1, tree2))


print(inputTree2, rf_dist, eucl_dist)



