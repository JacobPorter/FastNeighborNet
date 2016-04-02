package nnet;

//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.LinkedList;
//import java.util.List;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
//import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
//import java.util.concurrent.Future;
//import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Implements neighbor net
 * @version $Id:
 *
 * @author David Bryant
 * Adapted to Java by Daniel Huson 1.03
 * Modified by Jacob Porter.  Integrated into FastNN by Jacob Porter.
 *
 */
/**
 * Implements Neighbor Net method of Bryant and Moulton (2004).
 */
public class NeighborNetCanonical extends NetMakerOriginal {
	
	private volatile RowMinimum bestRM = null;
	//private final AtomicInteger[] myCounter = new AtomicInteger[] {new AtomicInteger(0), new AtomicInteger(0)};

    public NeighborNetCanonical(double[][] d, int numTaxa, int numThreads, ExecutorService pool) {
    	super(d, numTaxa, numThreads, pool);
    }
    
    private final int[] getIndices(int x) {
    	int row = (int) Math.ceil((-1 + Math.sqrt(1 + 8*x))/2);
    	int col = (int) (x -1 - (row-1)*(row)/2.0);
    	return new int[] {row, col};
    }
    
    private final synchronized void updateBest(NetNode myCx, NetNode myCy, double myBest) {
    	if (myBest < bestRM.value) {
    		bestRM = new RowMinimum(myCx, myCy, myBest);
    	}
    }
    
    public class FindNodesMulti implements Runnable {
    	
    	final NetNode[] netNodes;
    	final int num_clusters;
    	final int begin;
    	final int end;
    	final int id;
    	FindNodesMulti(final NetNode[] netNodes, final int num_clusters, final int begin, final int end, final int id) {
    		this.netNodes = netNodes;
    		this.num_clusters = num_clusters;
    		this.begin = begin;
    		this.end = end;
    		this.id = id;
    	}
    	
    	
    	public final void run() {
    		double Dpq;
            double Qpq;
            NetNode myCx;
            NetNode myCy;
            myCx = myCy = null;
            double myBest = Double.MAX_VALUE;
    		int[] beginIndices = getIndices(begin);
    		int[] endIndices = getIndices(end);
    		int i_start = beginIndices[0];
    		int j_start = beginIndices[1];
    		final int i_end = endIndices[0];
    		final int j_end = endIndices[1];
    		final int i = i_start;
    		//System.out.println(i_start + " " + j_start + " " + i_end + " " + j_end);
    		for (; i_start <= i_end; i_start++) {
                if (i_start != i) {
                	j_start = 0;
                }
            	if (i_start >= i_end && j_start >= j_end) {
            		//System.out.println("Breaking " + id);
            		break;
            	}
    			NetNode p = netNodes[i_start];
                if ((p.nbr != null) && (p.nbr.id < p.id)) /* We only evaluate one node per cluster */
                    continue;

                for (; j_start != i_start; j_start++) {
                	//System.out.println("inside j...");
                	if (i_start >= i_end && j_start >= j_end) {
                		//System.out.println("Breaking " + id);
                		break;
                	}
                	//myCounter[id].getAndIncrement();
                	NetNode q = netNodes[j_start];
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
                    if ((myCx == null || (Qpq < myBest)) && (p.nbr != q)) {
                        myCx = p;
                        myCy = q;
                        myBest = Qpq;
                    }
                }
            }
    		/* Update my best nodes */
    		//System.out.println("Ending values: " + i_start + " " + j_start + " " + id);
    		//System.out.println(bestRM.value);
    		if (myBest < bestRM.value) {
    			//System.out.println(myBest);
    			updateBest(myCx, myCy, myBest);
    		}
    	}
    }
    
//    private final int[] getThreadWork(int num_active) {
//    	
//    	return new int[] {0, 0};
//    }

    public final void findNodes(Stack<NetNode> amalgs, double D[][], NetNode[] netNodes, int num_nodes, int num_active, int num_clusters) {
        double Dpq;
        double Qpq;
    	//TODO: DB: 2x speedup by using symmetry
    	//TODO: JSP: Row sum adjustments can be done in linear time after it is initialized.  Done.
          Cx = Cy = null;
          /* Now minimize (m-2) D[C_i,C_k] - Sx - Sy */
          best = Double.MAX_VALUE;
          //myCounter[0].set(0);
         // myCounter[1].set(0);
          //System.out.println(num_active);
          if (numThreads == 1) {
	            for (int i = 0; i < num_active; i++) {
	            	NetNode p = netNodes[i];
	                if ((p.nbr != null) && (p.nbr.id < p.id)) /* We only evaluate one node per cluster */
	                    continue;
	                for (int j = 0; j != i; j++) {
	                	NetNode q = netNodes[j];
	                	//myCounter[0].getAndIncrement();
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
          } else {
        	  bestRM = new RowMinimum(null, null, Double.MAX_VALUE);
        	  final int amountOfWork = (num_active * (num_active - 1) / 2);
        	  final int workPerThread = amountOfWork / numThreads; // Math.floorDiv(amountOfWork, numThreads);
        	  int remainder = amountOfWork % numThreads;
        	  //System.out.println("work per thread and remainder: " + workPerThread + " " + remainder);
        	  List<Future> futures = new ArrayList<Future>(numThreads);
        	  int start = 1;
        	  for (int i = 0; i < numThreads; i++) {
        		  int end = start + workPerThread;
        		  if (remainder > 0) {
        			  end++;
        			  remainder--;
        		  }
        		  //Add stuff to queue.
        		  futures.add(pool.submit(new FindNodesMulti(netNodes,num_clusters, start, end, i)));
        		  start = end;
        	  }
        	  for(Future future : futures) {
					try {
					   future.get();
					} catch (Exception e) {}
				}
          	Cx = bestRM.me;
          	Cy = bestRM.row;
          	best = bestRM.value;
          }
          //System.out.println("Amount of work done this iteration: " + myCounter[0].get() + " " + myCounter[1].get());
          if (debug) {
        	  System.out.println("best: " + best);
        	  System.out.println(Cx);
        	  System.out.println(Cy);
        	  System.out.print("");
          }
          //return new NetNode[] {Cx, Cy};
    }
}
