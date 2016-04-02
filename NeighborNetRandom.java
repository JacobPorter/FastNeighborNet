package nnet;
import java.util.ArrayList;
import java.util.Stack;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;


public class NeighborNetRandom extends NetMakerOriginal {
	
	protected static enum RandomAmount {
	     LOGN, N, NLOGN
	}
	
	final RandomAmount myAmount;
	final Random myRandom;
	final int multiplier;
	int iterCount;
	volatile double myBest;
	volatile NetNode myCx;
	volatile NetNode myCy;
	
    public NeighborNetRandom(double[][] d, int numTaxa, int numThreads, ExecutorService pool, RandomAmount myAmount, int multiplier) {
    	super(d, numTaxa, numThreads, pool);
    	this.myAmount = myAmount;
		myRandom = ThreadLocalRandom.current();
		this.multiplier = multiplier;
    }
    
    private final int findSearchAmount(int total) {
		int amountToSearch;
		switch (myAmount) {
        case LOGN:
        	amountToSearch = (int) Math.ceil(Math.log10(total));
        	break;
        case N:
        	amountToSearch = total;
        	break;
        case NLOGN:
        	amountToSearch = (int) Math.ceil(Math.log10(total)) * total;
        	break;
        default:
        	amountToSearch = (int) Math.ceil(Math.log10(total)) * total;
        	break;
		}
		return multiplier * amountToSearch;
	}
    
    private synchronized void updateBest(NetNode myCx, NetNode myCy, double value) {
    	if (value < myBest) {
    		this.myCx = myCx;
    		this.myCy = myCy;
    		myBest = value;
    	}
    }
    
    class FindNodesMulti implements Runnable {
    	
		final int start;
		final int end;
		final NetNode[] netNodes;
		final int num_clusters;
		final int num_active;
    	FindNodesMulti(final NetNode[] netNodes, final int num_clusters, final int num_active, final int start, final int end) {
			this.start = start;
			this.end = end;
			this.netNodes = netNodes;
			this.num_clusters = num_clusters;
			this.num_active = num_active;
    	}
    	
    	public final void run() {
    		ThreadLocalRandom localRandom = ThreadLocalRandom.current();
    		int i = localRandom.nextInt(num_active);
	        int j;
	    	NetNode p, q;
	    	double Dpq;
	    	double Qpq;
	    	double best = Double.MAX_VALUE;
	    	NetNode Cx;
	    	NetNode Cy;
	    	Cx = Cy = null;
	    	for (int k = start; k < end; k++) {
	    		if (netNodes[i].nbr != null) {
	    			int iNbr = netNodes[i].nbr.positionID;
	    			j = localRandom.nextInt(num_active-2);
	    			if (i == j && num_active-1 == iNbr) {
	    				j = num_active-2;
	    			} else if (i == j && num_active-1 != iNbr) {
	    				j = num_active-1;
	    			} else if (iNbr == j && num_active-2 == i) {
	    				j = num_active-1;
	    			} else if (iNbr == j && num_active-2 != i) {
	    				j = num_active-2;
	    			}
	    		} else {
	    			j = localRandom.nextInt(num_active-1);
	    			if (i == j) {
		    			j = num_active-1;
		    		}
	    		}
	    		p = netNodes[i];
	    		q = netNodes[j];
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
	            i = j;
	    	}
	    	if (best < myBest) {
	    		updateBest(Cx, Cy, best);
	    	}
    	}
    }
    
    
    protected final void findNodes(Stack<NetNode> amalgs, double D[][], NetNode[] netNodes, int num_nodes, int num_active, int num_clusters) {
        best = Double.MAX_VALUE;
        final int searchAmount = findSearchAmount(num_active);
        if (numThreads == 1) {
	        int i = myRandom.nextInt(num_active);
	        int j;
	    	NetNode p, q;
	    	double Dpq;
	    	double Qpq;
	    	Cx = Cy = null;
	    	for (int k = 0; k < searchAmount ; k++) {
	    		if (netNodes[i].nbr != null) {
	    			int iNbr = netNodes[i].nbr.positionID;
	    			j = myRandom.nextInt(num_active-2);
	    			if (i == j && num_active-1 == iNbr) {
	    				j = num_active-2;
	    			} else if (i == j && num_active-1 != iNbr) {
	    				j = num_active-1;
	    			} else if (iNbr == j && num_active-2 == i) {
	    				j = num_active-1;
	    			} else if (iNbr == j && num_active-2 != i) {
	    				j = num_active-2;
	    			}
	    		} else {
	    			j = myRandom.nextInt(num_active-1);
	    			if (i == j) {
		    			j = num_active-1;
		    		}
	    		}
	    		p = netNodes[i];
	    		q = netNodes[j];
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
	            i = j;
	    	}
        } else {
			int start = 0;
			int work = searchAmount / numThreads; //Math.floorDiv(searchAmount, numThreads);
			int remainder = searchAmount % numThreads;
			myBest = Double.MAX_VALUE;
			ArrayList<Future<Integer>> futures = new ArrayList<Future<Integer>>();
			for (int i = 0; i < numThreads; i++) {
				int end = start + work + 1;
				if (remainder != 0) {
					end++;
					remainder--;
				}
				futures.add(pool.submit(new FindNodesMulti(netNodes, num_clusters, num_active, start, end), 0));
				start = end;
			}
			for(Future<Integer> future : futures) {
				try {
					future.get();
					} catch (Exception e) {}
				}
			Cx = myCx;
			Cy = myCy;
			best = myBest;
        }
    	//System.out.println(iterCount);
    	//System.out.println(num_active);
    	//System.out.println(Cx);
    	//System.out.println(Cy);
    }
    
}