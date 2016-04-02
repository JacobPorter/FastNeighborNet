package nnet;

import java.util.Stack;
import java.util.concurrent.ExecutorService;

import nnet.NetNode;

/**
 * Implements neighbor net
 * @version $Id:
 *
 * @author David Bryant
 * Adapted to Java by Daniel Huson 1.03
 * Modified by Jacob Porter.  Integrated into FastNN by Jacob Porter.
 *
 */
public abstract class NetMakerOriginal {
	
	enum NMMode {
		CANONICAL, RELAXED, RANDOM_N, RANDOM_NLOGN, RANDOM_LOGN, ORIGINAL //,FILTER
	}
	
    private double optionThreshold = 0.000001; // min weight of split that we consider
    private boolean makeSplits = true;
    private String optionVarianceName = "Ordinary_Least_Squares";
    private boolean optionConstrain = true;
    private int[] ordering = null; // the computed ordering
    public final static String DESCRIPTION = "Computes the Neighbor-Net network (Bryant and Moulton 2004)";
    protected final int ntax;
    protected final double[][] D;
    protected final int numThreads;
    protected final ExecutorService pool;
    volatile double best;
    volatile NetNode Cx;
    volatile NetNode Cy;
    //protected long count = 0;
    public final boolean debug = false;
    
	protected class RowMinimum {
		final NetNode me;
		final NetNode row;
		final double value;
		RowMinimum(final NetNode me, final NetNode row, final double value) {
			this.me = me;
			this.row = row;
			this.value = value;
		}
	}
    
    
    public NetMakerOriginal(double[][] d, int numTaxa, int numThreads, ExecutorService pool) {
    	ntax = numTaxa;
    	D = d;
    	this.numThreads = numThreads;
    	this.pool = pool;
    }

    /**
     * gets a short description of the algorithm
     *
     * @return a description
     */
    public String getDescription() {
        return DESCRIPTION;
    }

    /**
     * gets a cyclic ordering computed by the algorithm
     *
     * @return a cyclic ordering
     */
    public int[] getOrdering() {
        return ordering;
    }


    /**
     * Sets the power for least squares
     *
     * @param varName 0, 1 or 2
     */
    public void setOptionVariance(String varName) {
        this.optionVarianceName = varName;
    }

    /**
     * Gets the power for least squares
     *
     * @return the power
     */
    public String getOptionVariance() {
        return optionVarianceName;
    }

    public String selectVariance(String varianceName) {
        if (varianceName.equalsIgnoreCase("OrdinaryLeastSquares"))
            return "ols";
        else if (varianceName.equalsIgnoreCase("FitchMargoliash1"))
            return "fm1";
        else if (varianceName.equalsIgnoreCase("FitchMargoliash2"))
            return "fm2";
        else if (varianceName.equalsIgnoreCase("Estimated"))
            return "user";
        else
            return "ols"; //In case of uncertainty, do OLS
    }

    /**
     * Sets the constrained option for least squares
     *
     * @param flag set the constrained option?
     */
    public void setConstrain(boolean flag) {
        this.optionConstrain = flag;
    }

    /**
     * Gets the constrained option for least squares
     *
     * @return true, if will use the constrained least squares
     */
    public boolean getConstrain() {
        return optionConstrain;
    }
	
    /**
     * Run the neighbor net algorithm
     */
    public int[] runNeighborNet() { //(int ntax, double[][] D)  {


        //Special cases. When ntax<=3, the default circular ordering will work.
        if (ntax <= 3) // nnet can't do small data sets, so let's use split decomp
        {
            int[] ordering = new int[ntax + 1];
            for (int i = 0; i <= ntax; i++)
                ordering[i] = i;
            return ordering;

        }
        NetNode[] netNodes = new NetNode[ntax];

        /* Nodes are stored in a doubly linked list that we set up here */
        for (int i = ntax; i >= 1; i--) /* Initially, all singleton nodes are active */ {
            NetNode taxNode = new NetNode();
            taxNode.id = i;
            taxNode.positionID = i-1;
            netNodes[i-1] = taxNode;
            taxNode.distID = i-1;
        }
        
        
        /* Perform the agglomeration step */
        Stack<NetNode> amalgs = new Stack<NetNode>();
        int num_nodes = ntax;
        /*Initialize the row sums */
        //TODO: JSP -- Check correctness of row sums
        // Need to create heaps for filtered approach
        initialize(D, netNodes, num_nodes, num_nodes, num_nodes);
        num_nodes = agglomNodes(/*doc,*/ amalgs, D, netNodes, num_nodes);
        return expandNodes(/*doc,*/ num_nodes, ntax, amalgs, netNodes);
    }
    
    protected void initialize(double D[][], NetNode[] netNodes, int num_nodes, int num_active, int num_clusters) {
    	for (NetNode p : netNodes) { //(NetNode p = netNodes.next; p != null; p = p.next) {
            //NetNode p = netNodes[i]; 
        	if (p.nbr == null || p.nbr.id > p.id) {
        		for (int j = p.positionID+1; j < num_nodes; j++) {
        			NetNode q = netNodes[j];
                    if (q.nbr == null || (q.nbr.id > q.id) && (q.nbr != p)) {
                        Double Dpq = 0.0;
                        if ((p.nbr == null) && (q.nbr == null))
                            Dpq = D[p.distID][q.distID];
                        else if ((p.nbr != null) && (q.nbr == null))
                            Dpq = (D[p.distID][q.distID] + D[p.nbr.distID][q.distID]) / 2.0;
                        else if ((p.nbr == null) && (q.nbr != null))
                            Dpq = (D[p.distID][q.distID] + D[p.distID][q.nbr.distID]) / 2.0;
                        else
                            Dpq = (D[p.distID][q.distID] + D[p.distID][q.nbr.distID] + D[p.nbr.distID][q.distID] + D[p.nbr.distID][q.nbr.distID]) / 4.0;

                        p.Sx += Dpq;
                        if (p.nbr != null)
                            p.nbr.Sx += Dpq;
                        q.Sx += Dpq;
                        if (q.nbr != null)
                            q.nbr.Sx += Dpq;
                    }
                }
            }
        }
    }
    
    abstract protected void findNodes(Stack<NetNode> amalgs, double D[][], NetNode[] netNodes, int num_nodes, int num_active, int num_clusters);
    //abstract protected void initialize();
    //abstract protected void cleanup();
    
    protected final void findNodesDefault(Stack<NetNode> amalgs, double D[][], NetNode[] netNodes, int num_nodes, int num_active, int num_clusters) {
    	double Dpq;
        double Qpq;
    	//TODO: DB: 2x speedup by using symmetry
    	//TODO: JSP: Row sum adjustments can be done in linear time after it is initialized.  Done.

          Cx = Cy = null;
          /* Now minimize (m-2) D[C_i,C_k] - Sx - Sy */
          best = Double.MAX_VALUE;
          
          //System.out.println(num_active);
	            for (int i = 0; i < num_active; i++) {
	            	NetNode p = netNodes[i];
	                if ((p.nbr != null) && (p.nbr.id < p.id)) /* We only evaluate one node per cluster */
	                    continue;
	                for (int j = 0; j != i; j++) {
	                	NetNode q = netNodes[j];
	                    if ((q.nbr != null) && (q.nbr.id < q.id)) /* We only evaluate one node per cluster */
	                        continue;
	                    if (q.nbr == p) /* We only evaluate nodes in different clusters */
	                        continue;
	                    if ((p.nbr == null) && (q.nbr == null))
	                        Dpq = D[p.distID][q.distID];
	                    else if ((p.nbr != null) && (q.nbr == null))
	                        Dpq = (D[p.distID][q.distID] + D[p.nbr.distID][q.distID]) / 2.0;
	                    else if ((p.nbr == null) && (q.nbr != null))
	                        Dpq = (D[p.distID][q.distID] + D[p.distID][q.nbr.distID]) / 2.0;
	                    else
	                        Dpq = (D[p.distID][q.distID] + D[p.distID][q.nbr.distID] + D[p.nbr.distID][q.distID] + D[p.nbr.distID][q.nbr.distID]) / 4.0;
	                    Qpq = ((double) num_clusters - 2.0) * Dpq - p.Sx - q.Sx;
	                    //count++;
	                    /* Check if this is the best so far */
	                    if ((Cx == null || (Qpq < best)) && (p.nbr != q)) {
	                        Cx = p;
	                        Cy = q;
	                        best = Qpq;
	                    }
	                }
	            }
    }
    
    /**
     * Expands the net nodes to obtain the ordering, quickly
     *
     * @param num_nodes number of nodes
     * @param ntax      number of taxa
     * @param amalgs    stack of amalagations
     * @param netNodes  the net nodes
     */
    private final int[] expandNodes(/*Document doc,*/ int num_nodes, int ntax, Stack<NetNode> amalgs, NetNode[] netNodes) { // throws CanceledException {

        int[] ordering = new int[ntax + 1];  //Does this need to be changed??
        //System.err.println("expandNodes");
        NetNode x, y, z, u, v, a;
        
        /* Set up the circular order for the first three nodes */
        //x = netNodes.next;
        //y = x.next;
        //z = y.next;
        x = netNodes[0];
        y = netNodes[1];
        z = netNodes[2];

        x.next = y;
        y.next = z;
        z.next = x;
        x.prev = z;
        y.prev = x;
        z.prev = y;
        if (debug) {
        	System.out.println("expand Nodes.  Last nodes");
        	System.out.println(x);
        	System.out.println(y);
        	System.out.println(z);
        }

        /* Now do the rest of the expansions */
        while (!amalgs.empty()) {
        	/* Find the three elements replacing u and v. Swap u and v around if v comes before u in the
          circular ordering being built up */
            u = (NetNode) (amalgs.pop());
            // System.err.println("POP: u="+u);
            v = u.nbr;
            x = u.ch1;
            y = u.ch2;
            z = v.ch2;
            if (v != u.next) {
                NetNode tmp = u;
                u = v;
                v = tmp;
                tmp = x;
                x = z;
                z = tmp;
            }

            /* Insert x,y,z into the circular order */
            x.prev = u.prev;
//            try {
            x.prev.next = x;
//            } catch (NullPointerException e) {
//            	System.out.println("Null pointer found...");
//            	System.out.println(x);
//            	System.out.println(amalgs.size());
//            	throw e;
//            }
            x.next = y;
            y.prev = x;
            y.next = z;
            z.prev = y;
            z.next = v.next;
            z.next.prev = z;
        }

        /* When we exit, we know that the point x points to a node in the circular order */
        /* We loop through until we find the node after taxa zero */
        while (x.id != 1) {
            x = x.next;
        }

        /* extract the ordering */
        a = x;
        int t = 0;
        do {
            // System.err.println("a="+a);
            ordering[++t] = a.id;
            a = a.next;
        } while (a != x);
        return ordering;
    }
    

    /**
     * Agglomerates the nodes
     */
    protected int agglomNodes(/*Document doc,*/ Stack<NetNode> amalgs, double D[][], NetNode[] netNodes, int num_nodes) { //throws CanceledException {
        //System.err.println("agglomNodes");
        //double best;
        int num_active = num_nodes;
        int num_clusters = num_nodes;
        //final Object lock1 = new Object();
        //best = Double.MAX_VALUE;
        
        while (num_active > 3) {

            /* Special case
            If we let this one go then we get a divide by zero when computing Qpq */
            if (num_active == 4 && num_clusters == 2) {
            	NetNode q;
                //NetNode p = netNodes.next;
            	NetNode p = netNodes[0];
            	if (p.nbr != netNodes[1]) {
            		q = netNodes[1];
            	} else {
            		q = netNodes[2];
            	}
                if (D[p.distID][q.distID] + D[p.nbr.distID][q.nbr.distID] < D[p.distID][q.nbr.distID] + D[p.nbr.distID][q.distID]) {
                    agg3way(p, q, q.nbr, amalgs, D, netNodes, num_nodes, num_active);
                    num_nodes += 2;
                } else {
                    agg3way(p, q.nbr, q, amalgs, D, netNodes, num_nodes, num_active);
                    num_nodes += 2;
                }
                break;
            }
            if (num_active <= 1024) {
            	//System.err.println("Here!!");
        		findNodesDefault(amalgs, D, netNodes, num_nodes, num_active, num_clusters);
        	} else {
        		findNodes(amalgs, D, netNodes, num_nodes, num_active, num_clusters);
        	}
//            System.out.println("Cx: " + Cx);
//            if (Cx.nbr != null) {
//            	System.out.println("Cx.nbr: " + Cx.nbr);
//            }
//            System.out.println("Cy: " + Cy);
//            if (Cy.nbr != null) {
//            	System.out.println("Cy.nbr: " + Cy.nbr);
//            }
            //System.out.print(num_active + " " + num_clusters + "\n");
            if (Cx.id > Cy.id) {  //TODO: If the order is different, the row sums will be different.  Bug?
            	NetNode temp = Cx;
            	Cx = Cy;
            	Cy = temp;
            }
//            System.err.println("Best pair found: " + best + " " + Cx + " " + Cy);
//            for (NetNode nn : netNodes) {
//            	if (nn == null) {
//            		break;
//            	}
//            	System.err.print(nn.id + " ");
//            }
//            System.err.print("\n");
            int[] updates = handleAgglomerationEvent(Cx, Cy, amalgs, D, netNodes, num_nodes, num_active, num_clusters);
            num_nodes = updates[0];
            num_active = updates[1];
            num_clusters = updates[2];
        }
        return num_nodes;
    }
    
    protected final int[] handleAgglomerationEvent(NetNode Cx, NetNode Cy, Stack<NetNode> amalgs, 
    		double D[][], NetNode[] netNodes, int num_nodes, int num_active, int num_clusters) {
    	NetNode x;
    	NetNode y;
        int m;
        Double Qpq;
    	/* Find the node in each cluster */
        x = Cx;
        y = Cy;
        //System.out.println(x);
        //System.out.println(y);
        /*JSP row sums in the node can be replaced with variables since they are not needed. */
        Double Cx_Rx=0.0;
        Double Cx_nbr_Rx=0.0;
        Double Cy_Rx=0.0;
        Double Cy_nbr_Rx=0.0;
        if (Cx.nbr != null || Cy.nbr != null) {
            Cx_Rx = ComputeRx(Cx, Cx, Cy, D, netNodes, num_active);
            if (Cx.nbr != null)
                Cx_nbr_Rx = ComputeRx(Cx.nbr, Cx, Cy, D, netNodes, num_active);
            Cy_Rx = ComputeRx(Cy, Cx, Cy, D, netNodes, num_active);
            if (Cy.nbr != null)
                Cy_nbr_Rx = ComputeRx(Cy.nbr, Cx, Cy, D, netNodes, num_active);
        }

        m = num_clusters;
        if (Cx.nbr != null)
            m++;
        if (Cy.nbr != null)
            m++;

        best = ((double) m - 2.0) * D[Cx.distID][Cy.distID] - Cx_Rx - Cy_Rx;
        if (Cx.nbr != null) {
            Qpq = ((double) m - 2.0) * D[Cx.nbr.distID][Cy.distID] - Cx_nbr_Rx - Cy_Rx;
            if (Qpq < best) {
                x = Cx.nbr;
                y = Cy;
                best = Qpq;
            }
        }
        if (Cy.nbr != null) {
            Qpq = ((double) m - 2.0) * D[Cx.distID][Cy.nbr.distID] - Cx_Rx - Cy_nbr_Rx;
            if (Qpq < best) {
                x = Cx;
                y = Cy.nbr;
                best = Qpq;
            }
        }
        if ((Cx.nbr != null) && (Cy.nbr != null)) {
            Qpq = ((double) m - 2.0) * D[Cx.nbr.distID][Cy.nbr.distID] - Cx_nbr_Rx - Cy_nbr_Rx;
            if (Qpq < best) {
                x = Cx.nbr;
                y = Cy.nbr;
                best = Qpq;
            }
        }
        /* Subtract old cluster distances */
        NetNode u;
        for (int i = 0; i < num_active; i++) {
        	NetNode p = netNodes[i];
        	if (i != x.positionID && i != y.positionID) {
        		subtractClusterDistance(p, x);
        		subtractClusterDistance(p, y);
        	}
        }
        if ((null == x.nbr) && (null == y.nbr)) {   /* Both vertices are isolated...add edge {x,y} */
            u = agg2way(x, y);
            num_clusters--;
        } else if (null == x.nbr) {     /* X is isolated,  Y  is not isolated*/
            u = agg3way(x, y, y.nbr, amalgs, D, netNodes, num_nodes, num_active);
            num_nodes += 2;
            num_active--;
            num_clusters--;
            x.positionID = -1;
            y.positionID = -1;
            y.nbr.positionID = -1;
            if (y.nbr != null) {y.nbr.positionID = -1;}
        } else if ((null == y.nbr) || (num_active == 4)) { /* Y is isolated,  X is not isolated
                                                    OR theres only four active nodes and none are isolated */
            u = agg3way(y, x, x.nbr, amalgs, D, netNodes, num_nodes, num_active);
            num_nodes += 2;
            num_active--;
            num_clusters--;
            x.positionID = -1;
            y.positionID = -1;
            x.nbr.positionID = -1;
        } else {  /* Both nodes are connected to others and there are more than 4 active nodes */
            u = agg4way(x.nbr, x, y, y.nbr, amalgs, D, netNodes, num_nodes, num_active);
            num_nodes += 4;
            num_active -= 2;
            num_clusters--;
        }
        /* Add new cluster distances */
        updateClusterDistances(u, D, netNodes, num_active);
//        u.Sx = 0;
//        u.nbr.Sx = 0;
//        for (int i = 0; i < num_active; i++) {
//        	NetNode p = netNodes[i];
//        	if  ((p.nbr == null || p.nbr.id > p.id) && (u.nbr != p) && (u != p)) {
//        		Double Dpu = 0.0;
//                if (p.nbr == null) {// && (u.nbr != null))
//                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID]) / 2.0;
//                } else {
//                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID] + D[p.nbr.distID][u.distID] + D[p.nbr.distID][u.nbr.distID]) / 4.0;
//                }
//                // add Dpu to a hashmap for the filtered approach
//                p.Sx += Dpu;
//                if (p.nbr != null) {p.nbr.Sx += Dpu; }  // Could remove addition and set to p.Sx
//                u.Sx += Dpu;
//        	}
//        }
//        u.nbr.Sx = u.Sx;
        // node u needs to be added to the heaps...
        if (debug) {
        	System.out.println("New node=" + u);
        	System.out.print("");
        }
    	return new int[] {num_nodes, num_active, num_clusters};
    }
    
    protected void updateClusterDistances(NetNode u, double D[][], NetNode[] netNodes, int num_active) {
    	u.Sx = 0;
        u.nbr.Sx = 0;
        for (int i = 0; i < num_active; i++) {
        	NetNode p = netNodes[i];
        	if  ((p.nbr == null || p.nbr.id > p.id) && (u.nbr != p) && (u != p)) {
        		Double Dpu = 0.0;
                if (p.nbr == null) {// && (u.nbr != null))
                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID]) / 2.0;
                } else {
                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID] + D[p.nbr.distID][u.distID] + D[p.nbr.distID][u.nbr.distID]) / 4.0;
                }
                // add Dpu to a hashmap for the filtered approach
                p.Sx += Dpu;
                if (p.nbr != null) {p.nbr.Sx += Dpu; }  // Could remove addition and set to p.Sx
                u.Sx += Dpu;
        	}
        }
        u.nbr.Sx = u.Sx;
    }
    

    /**
     * Computes the Rx
     *
     * @param z        a node
     * @param Cx       a node
     * @param Cy       a node
     * @param D        the distances
     * @param netNodes the net nodes
     * @return the Rx value
     */
    final protected double ComputeRx(NetNode z, NetNode Cx, NetNode Cy, double[][] D,
                                    NetNode[] netNodes, int num_active) {
        double Rx = 0.0;

        for (int i = 0; i < num_active; i++) {
        	NetNode p = netNodes[i];
            if (p == Cx || p == Cx.nbr || p == Cy || p == Cy.nbr || p.nbr == null)
                Rx += D[z.distID][p.distID];
            else /* p.nbr != null */
                Rx += D[z.distID][p.distID] / 2.0; /* We take the average of the distances */ //JSP: Why??
        }
        return Rx;
    }
    

    /**
     * agglomerate 2 nodes
     *
     * @param x one node
     * @param y other node
     */
    final protected NetNode agg2way(NetNode x, NetNode y) {
    	if (debug) {
    		System.out.println("2way");
    	}
        x.nbr = y;
        y.nbr = x;
        return x;
    }

    /**
     * agglomerate 3 nodes.
     * Note that this version doesn't update num_nodes, you need to
     * num_nodes+=2 after calling this!
     * x and z are always in different clusters
     * @param x one node
     * @param y other node
     * @param z other node
     * @return one of the new nodes
     */
    final protected NetNode agg3way(NetNode x, NetNode y, NetNode z,
                                   Stack<NetNode> amalgs, double[][] D, NetNode[] netNodes, int num_nodes, int num_active) {
/* Agglomerate x,y, and z to give TWO new nodes, u and v */
/* In terms of the linked list: we replace x and z
  	 by u and v and remove y from the linked list.
  	 and replace y with the new node z
    Returns a pointer to the node u */
//printf("Three way: %d, %d, and %d\n",x.id,y.id,z.id);
    	
    	/* Linearly update cluster distances by subtracting off cluster distances to x and z. */
    	//TODO: Could there be a subtle bug in the update?
//    	for (int i = 0; i < num_active; i++) {
//        	NetNode p = netNodes[i];
//        	subtractClusterDistance(p, x);
//        	subtractClusterDistance(p, z);
//        }
    	if (debug) {
    		System.out.println("3way");
    	}
        NetNode u = new NetNode();
        u.id = num_nodes + 1;
        u.ch1 = x;
        u.ch2 = y;

        NetNode v = new NetNode();
        v.id = num_nodes + 2;
        v.ch1 = y;
        v.ch2 = z;

        // JSP: Need to invalidate x, y, z by setting positionID to -1
        /* Replace x by u in the linked list */
        if (x.positionID == -1) {
        	System.err.println(x);
        }
        netNodes[x.positionID] = u;
        u.positionID = x.positionID;
        u.distID = x.distID;
        // Needed for the filtered approach
        
        
        /* Replace z by v in the linked list */
        netNodes[z.positionID] = v;
        v.positionID = z.positionID;
        v.distID = z.distID;
        // Needed for the filtered approach
        
        

        /* Remove y from the linked list */
        if (y.positionID == -1) {
        	System.err.println(y);
        }
        netNodes[y.positionID] = netNodes[num_active-1];
        netNodes[y.positionID].positionID = y.positionID;
        netNodes[num_active-1] = null;
        // Needed for the filtered approach

        /* Add an edge between u and v, and add u into the list of amalgamations */
        u.nbr = v;
        v.nbr = u;
        

        /* Update distance matrix */
        /* Linearly update cluster row sums */
        for (int i = 0; i < num_active-1; i++) {  //Subtract off the old value.  Add the new value.
        	NetNode p = netNodes[i];  //What about the 4-way? //What if p = u or v?
        	D[u.distID][p.distID] = D[p.distID][u.distID] = (2.0 / 3.0) * D[x.distID][p.distID] + D[y.distID][p.distID] / 3.0;
            D[v.distID][p.distID] = D[p.distID][v.distID] = (2.0 / 3.0) * D[z.distID][p.distID] + D[y.distID][p.distID] / 3.0;
//            if ((p.nbr == null || p.nbr.id > p.id) && (u.nbr != p) && (u != p)) {
//                Double Dpu = 0.0;
//                if (p.nbr == null) {// && (u.nbr != null))
//                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID]) / 2.0;
//                } else {
//                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID] + D[p.nbr.distID][u.distID] + D[p.nbr.distID][u.nbr.distID]) / 4.0;
//                }
//                p.Sx += Dpu;
//                if (p.nbr != null) {p.nbr.Sx += Dpu; }  // Could remove addition and set to p.Sx
//                u.Sx += Dpu;
//            }
        }
        // What about the distance between u and v?
        D[u.distID][u.distID] = D[v.distID][v.distID] = 0.0;
        //v.Sx = u.Sx;
        amalgs.push(u);
        return u;
    }
    
    /**
     * Subtracts the cluster distance from p to x from p.Sx
     * @param p A NetNode to subtract from
     * @param x A NetNode to calcualte the cluster distance to. 
     */
    private final void subtractClusterDistance(NetNode p, NetNode x) {
    	if (p != x && p != x.nbr && (p.nbr == null || (p.nbr.id > p.id))) {
        	Double Dpx = 0.0;
            if ((p.nbr == null) && (x.nbr == null))
                Dpx = D[p.distID][x.distID];
            else if ((p.nbr != null) && (x.nbr == null))
                Dpx = (D[p.distID][x.distID] + D[p.nbr.distID][x.distID]) / 2.0;
            else if ((p.nbr == null) && (x.nbr != null))
                Dpx = (D[p.distID][x.distID] + D[p.distID][x.nbr.distID]) / 2.0;
            else
                Dpx = (D[p.distID][x.distID] + D[p.distID][x.nbr.distID] + D[p.nbr.distID][x.distID] + D[p.nbr.distID][x.nbr.distID]) / 4.0;
            p.Sx -= Dpx;
            if (p.nbr != null)
                p.nbr.Sx -= Dpx;
    	}
    }

    /**
     * Agglomerate four nodes
     *
     * @param x2 a node
     * @param x  a node
     * @param y  a node
     * @param y2 a node
     * @return the new number of nodes
     */
    final protected NetNode agg4way(NetNode x2, NetNode x, NetNode y, NetNode y2,
                               Stack<NetNode> amalgs, double[][] D, NetNode[] netNodes, int num_nodes, int num_active) {
    	/* Replace x2,x,y,y2 by with two vertices... performed using two
  	 	3 way amalgamations */
    	if (debug) {
    		System.out.println("4way");
    	}
    	NetNode u, v;
        u = agg3way(x2, x, y, amalgs, D, netNodes, num_nodes, num_active); /* Replace x2,x,y by two nodes, equal to x2_prev.next and y_prev.next. */
        num_nodes += 2;
        v = agg3way(u, u.nbr, y2, amalgs, D, netNodes, num_nodes, num_active-1); /* z = y_prev . next */
        num_nodes += 2;
        x2.positionID = -1;
        x.positionID = -1;
        y.positionID = -1;
        y2.positionID = -1;
        u.positionID = -1;
        u.nbr.positionID = -1;
        return v;
    }

}