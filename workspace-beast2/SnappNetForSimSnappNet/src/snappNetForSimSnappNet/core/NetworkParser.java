package snappNetForSimSnappNet.core;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.core.Input.Validate;

/**
 * Parse the network of extended Newick format.
 * @author Huw Ogilvie
 */

@Description("Parse the network of extended Newick format.")
public class NetworkParser extends Network implements StateNodeInitialiser {
    public final Input<Network> networkInput = new Input<>("initial", "Network to initialize.");
    public final Input<Tree> treeInput =
            new Input<>("tree", "Tree initialized from extended newick string.", Validate.REQUIRED);
    public final Input<Boolean> adjustTipHeightsInput =
            new Input<>("adjustTipHeights", "Whether tipHeights shall be adjusted (default is true).", true);

    private List<String> leafOrder;
    private int nextSpeciationNr;
    private int nextReticulationNr;

    public NetworkParser() {
    }

    public NetworkParser(final Tree tree) {
	Log.debug.println("Je construit mon Network Parser depuis un arbre");//CE
        treeInput.setValue(tree, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        final Tree tree = treeInput.get();
        final Node treeRoot = tree.getRoot();

        // Step (1) is to initialize the node counts and array
        leafOrder = new ArrayList<>();
        leafNodeCount = 0;
        speciationNodeCount = 0;
        int hybridNodeCount = 0;	
	
	
        for (Node n: tree.getNodesAsArray()) {
            if (n.getID() != null && n.getID().startsWith("#H")) {
                hybridNodeCount++;
            } else if (n.isLeaf()) {
            	leafOrder.add(n.getID());
                leafNodeCount++;
            } else if (!n.isRoot()) {
                speciationNodeCount++;
            }
        }

        leafOrder.sort(null);

        assert hybridNodeCount % 2 == 0;
        reticulationNodeCount = hybridNodeCount / 2;
        nodeCount = leafNodeCount + speciationNodeCount + reticulationNodeCount + 1;
        nodes = new NetworkNode[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            nodes[i] = new NetworkNode(this); 
        }

        nextSpeciationNr = leafNodeCount;
        nextReticulationNr = leafNodeCount + speciationNodeCount;
	 
        // Step (2) is to recursively copy the tree to the network
        rebuildNetwork(treeRoot);
	 
        // Update the cached parents and children for each node
        updateRelationships();

        // Step (3) adjust network tip height to ZERO
        if (adjustTipHeightsInput.get()) {
            // all nodes should be at zero height if no date-trait is available
            for (NetworkNode tip: getLeafNodes()) {
                tip.setHeight(0.0);
            }
        }

        super.initAndValidate();       
        
	 
    }

    private Integer rebuildNetwork(final Node treeNode) {

	    Integer branchNumber;
        NetworkNode newNode;

        final String nodeLabel = treeNode.getID();
        final double nodeHeight = treeNode.getHeight();
        final int matchingNodeNr = getNodeNumber(nodeLabel);		
	 
        if (matchingNodeNr < 0) {
            int newNodeNumber;
            double inheritProb = 0.5;

            if (treeNode.isRoot()) {
                newNodeNumber = nodeCount - 1;
            } else if (nodeLabel != null && nodeLabel.startsWith("#H")) {
                if (treeNode.getMetaDataNames().contains("gamma"))
                    inheritProb = (Double) treeNode.getMetaData("gamma");
                newNodeNumber = nextReticulationNr;
                nextReticulationNr++;
            } else if (treeNode.isLeaf()) {
                newNodeNumber = leafOrder.indexOf(nodeLabel);
            } else {
                newNodeNumber = nextSpeciationNr;
                nextSpeciationNr++;
            }

            newNode = nodes[newNodeNumber];
            newNode.label = nodeLabel;
            newNode.height = nodeHeight;
            newNode.inheritProb = inheritProb;

            branchNumber = getBranchNumber(newNodeNumber);
	 
        } else {	   
            newNode = nodes[matchingNodeNr];
            if (treeNode.getMetaDataNames().contains("gamma")){ 
                newNode.inheritProb = 1.0 - (double)treeNode.getMetaData("gamma");
	    }

            branchNumber = getBranchNumber(matchingNodeNr) + 1;
        }

        for (Node child: treeNode.getChildren()) {
            final Integer childBranchNr = rebuildNetwork(child);
            newNode.childBranchNumbers.add(childBranchNr);
        }

        return branchNumber;
    }

    @Override
    public void initStateNodes() {
        if (networkInput.get() != null) {
            networkInput.get().assignFrom(this);
        }
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
        if (networkInput.get() != null) {
            stateNodes.add(networkInput.get());
        }
    }
}
