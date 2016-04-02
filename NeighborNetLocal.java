package nnet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.concurrent.ThreadLocalRandom;


class NeighborNetLocal extends NetMakerOriginal {
	
	final Random myRandom;
	final int[] rowPermutation;
	boolean firstTime = true;
	int top;
	final boolean additive;
	//volatile double myBest;
	//final ArrayList<Integer> counter = new ArrayList<Integer>(5*ntax);
	//int count = 0;
	//final ArrayList<Integer> rowPermutation = new ArrayList<Integer>(ntax);
	
	public NeighborNetLocal(double[][] d, int numTaxa, int numThreads, boolean additive, ExecutorService pool) {
		super(d, numTaxa, numThreads, pool);
		//myRandom = new Random(System.currentTimeMillis());
		this.additive = additive;
		myRandom = ThreadLocalRandom.current();
		rowPermutation = new int[numTaxa];
	}
	
	
	class FindRowMinMulti implements Callable<ArrayList<RowMinimum>> {
		
		final int start;
		final int end;
		final NetNode[] netNodes;
		final int num_clusters;
		final NetNode p;
		FindRowMinMulti(final NetNode[] netNodes, final int num_clusters, final NetNode p, final int start, final int end) {
			this.start = start;
			this.end = end;
			this.netNodes = netNodes;
			this.num_clusters = num_clusters;
			this.p = p;
		}
		
		final public ArrayList<RowMinimum> call() {
			double Dpq;
			double Qpq;
			double myMin = Double.MAX_VALUE;
			ArrayList<RowMinimum> myMinimums = new ArrayList<RowMinimum>(2);
			for (int row = start; row < end; row++) {
				//NetNode p = netNodes[me];
				NetNode q = netNodes[row];
				if ((p == q) || ((p.nbr != null) && (p.nbr == q))) {
					continue;
				}
				//calculate distance
				//Parallelization: put work on a work queue?
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
				if (Qpq < myMin) {
					myMin = Qpq;
					myMinimums.clear();
					myMinimums.add(new RowMinimum(p, q, Qpq));
				} else if (Qpq == myMin) {
					myMinimums.add(new RowMinimum(p, q, Qpq));
				}
			}
			return myMinimums;
		}
	}
	
	

	
	final private List<RowMinimum> findRowMin(final NetNode p, HashMap<NetNode, List<RowMinimum>> foundRowMinimums, NetNode[] netNodes, int num_active, int num_clusters) {
		if (foundRowMinimums.containsKey(p)) {
			return foundRowMinimums.get(p);
		}
		if (p.nbr != null && foundRowMinimums.containsKey(p.nbr)) {
			return foundRowMinimums.get(p.nbr);
		}
		ArrayList<RowMinimum> myMinimums = new ArrayList<RowMinimum>(3);
		double myMin = Double.MAX_VALUE;
		if ((numThreads == 1 ) || (num_active <= 10000)) {
			double Dpq;
			double Qpq;
			for (int row = 0; row < num_active; row++) {
				//NetNode p = netNodes[me];
				NetNode q = netNodes[row];
				if ((p == q) || ((p.nbr != null) && (p.nbr == q))) {
					continue;
				}
				//calculate distance
				//Parallelization: put work on a work queue?
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
				if (Qpq < myMin) {
					myMin = Qpq;
					myMinimums.clear();
					myMinimums.add(new RowMinimum(p, q, Qpq));
				} else if (Qpq == myMin) {
					myMinimums.add(new RowMinimum(p, q, Qpq));
				}
				
			}
		} else {
			int start = 0;
			int work = num_active / numThreads; //Math.floorDiv(num_active, numThreads);
			int remainder = num_active % numThreads;
			ArrayList<Future<ArrayList<RowMinimum>>> futures = new ArrayList<Future<ArrayList<RowMinimum>>>();
			for (int i = 0; i < numThreads; i++) {
				int end = start + work + 1;
				if (remainder != 0) {
					end++;
					remainder--;
				}
				futures.add(pool.submit(new FindRowMinMulti(netNodes, num_clusters, p, start, end)));
				start = end;
			}
			for(Future<ArrayList<RowMinimum>> future : futures) {
				try {
						ArrayList<RowMinimum> al = future.get();
						if (al.get(0).value < myMin) {
							myMin = al.get(0).value;
							myMinimums = al;
						} else if (al.get(0).value == myMin) {
							for (RowMinimum rm : al) {
								myMinimums.add(rm);
							}
						}
					} catch (Exception e) {}
				}
			}
		foundRowMinimums.put(p, myMinimums);
		return myMinimums;
	}
	
	private final void swap(int[] x, int a, int b) {
		int t = x[a];
		x[a] = x[b];
		x[b] = t;
	}
	
	
	/**
	 *  Sampling without replacement based on Algorithm 3.4.2S of Knuth's book Seminumeric Algorithms.
	 *  http://stackoverflow.com/questions/311703/algorithm-for-sampling-without-replacement
	 */
	final protected void findNodes(final Stack<NetNode> amalgs, final double D[][], NetNode[] netNodes, final int num_nodes, final int num_active, final int num_clusters) {
		HashMap<NetNode, List<RowMinimum>> foundRowMinimums = new HashMap<NetNode, List<RowMinimum>>(); //use clear?
		ArrayList<RowMinimum> myMinimums = new ArrayList<RowMinimum>(5);
		RowMinimum combineMe = null;
		List<RowMinimum> testRowMin;
		List<RowMinimum> testOtherRow;
		if (firstTime) {
			for (int i = 0; i < ntax; i++) {
				rowPermutation[i] = i;
			}
			firstTime = false;
			top = ntax - 1;
		}
		//System.err.println("Performing swapping.  Num active: " + num_active + " top: " + top);
		outerloop:
		for (int i=top+1; i>0; i--) {
			int swapCell = myRandom.nextInt(i);
			if (rowPermutation[swapCell] >= num_active) {
				swap(rowPermutation, swapCell, top);
				if (i == top+1) {
					i--;
				} else {
					i++;
				}
				top--;
				//swap(rowPermutation, swapCell, i-1);
				continue;
			}
			swap(rowPermutation, i-1, swapCell);
			//System.err.print(rowPermutation[i-1] + " ");
			NetNode p = netNodes[rowPermutation[i-1]];
			if ((p.nbr != null) && (p.nbr.id < p.id)) { /* We only evaluate one node per cluster */
                continue;
			}
			testRowMin = findRowMin(p, foundRowMinimums, netNodes, num_active, num_clusters);
			for (RowMinimum myRM : testRowMin) {
				testOtherRow = findRowMin(myRM.row, foundRowMinimums, netNodes, num_active, num_clusters);
				for (RowMinimum testRM : testOtherRow) {
					//if (testRM.row.positionID == i) { //Change to test neighbor
					if ((testRM.row == p) ||  ((testRM.row.nbr != null) && (testRM.row.nbr == p)) || 
							((testRM.row.nbr != null) && (p.nbr != null) && (testRM.row.nbr == p.nbr)) ||
							((p.nbr != null) && (testRM.row == p.nbr)  )  )   { //Is all this necessary?
						myMinimums.add(testRM);
						break;  //Needed
					}
				}
			}
			//Check if local minimums have been found
			if (!myMinimums.isEmpty()) {
				int choice = myRandom.nextInt(myMinimums.size());
				combineMe = myMinimums.get(choice);
				Cx = combineMe.me;
				Cy = combineMe.row;
				if (additive) {
					int initial = choice;
					int max = myMinimums.size();
					while (true) {
						NetNode testNode = null;
						 for (int j = num_active-1; i > 0; i--) {
							 if (netNodes[j] == Cx || netNodes[j] == Cy || netNodes[j] == Cx.nbr || netNodes[j] == Cy.nbr) {
								 continue;
							 } else {
								 testNode = netNodes[j];
								 break;
							 }
						 }
						 //testNode needs to be the lexicographically least node among its neighbors
						 if (testNode.nbr != null && testNode.nbr.id < testNode.id) {
							 testNode = testNode.nbr;
						 }
						//Compute transformed distance between Cx and testNode
						 double originalDistCx = calculateDpq(Cx, testNode);
						 //double originalDistCy = calculateDpq(Cy, testNode);
						 double originalQ = ((double) num_clusters - 2.0) * originalDistCx - Cx.Sx - testNode.Sx;
						 //Compute transformed distance between Cx, Cy and testNode as if Cx and Cy were agglomerated.
						 double newQ = findAgglomeratedQ(Cx, Cy, testNode, num_clusters, netNodes, num_nodes, num_active);
						 //Compare distances.  If they are equal, then accept the join and agglomerate.  Otherwise DO NOT break the loop.
						 if (Math.abs(originalQ-newQ) < .0000001 ) {
							 break outerloop;
						 }
						 choice++;
						 choice %= max;
						 if (choice == initial) {
							 break;
						 }
					}
				} else {
					break outerloop;  //Needed
				}
			}
		}
		

		//System.err.print("\n");
	}
	
	private final double calculateDpq(NetNode p, NetNode q) {
		double Dpq;//, Qpq;
		if ((p.nbr == null) && (q.nbr == null))
            Dpq = D[p.distID][q.distID];
        else if ((p.nbr != null) && (q.nbr == null))
            Dpq = (D[p.distID][q.distID] + D[p.nbr.distID][q.distID]) / 2.0;
        else if ((p.nbr == null) && (q.nbr != null))
            Dpq = (D[p.distID][q.distID] + D[p.distID][q.nbr.distID]) / 2.0;
        else
            Dpq = (D[p.distID][q.distID] + D[p.distID][q.nbr.distID] + D[p.nbr.distID][q.distID] + D[p.nbr.distID][q.nbr.distID]) / 4.0;
        //Qpq = ((double) num_clusters - 2.0) * Dpq - p.Sx - q.Sx;
        return Dpq;
	}
	
	private final double findAgglomeratedQ(final NetNode Cx, final NetNode Cy, final NetNode testNode, 
			int num_clusters, NetNode[] netNodes, int num_nodes, int num_active) {
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
        double subtractedRowSum = testNode.Sx - calculateDpq(Cx, testNode) - calculateDpq(Cy, testNode);
        NetNode u;
//        for (int i = 0; i < num_active; i++) {
//        	NetNode p = netNodes[i];
//        	if (i != x.positionID && i != y.positionID) {
//        		subtractClusterDistance(p, x);
//        		subtractClusterDistance(p, y);
//        	}
//        }
        if ((null == x.nbr) && (null == y.nbr)) {   /* Both vertices are isolated...add edge {x,y} */
//            u = agg2way(x, y);
//            num_clusters--;
        	double Dpu = calculateClusterDistLocal(x, testNode);
        	double clusterRowSum = 0;
        	for (int i = 0; i < num_active; i++) {
        		NetNode p = netNodes[i];
        		if ((p.nbr == null || p.nbr.id > p.id) && (p != x && p != y)) {
        	        if (p.nbr == null) {// && (u.nbr != null))
        	        	clusterRowSum += (D[p.distID][x.distID] + D[p.distID][y.distID]) / 2.0;
        	        } else {
        	        	clusterRowSum += (D[p.distID][x.distID] + D[p.distID][y.distID] + D[p.nbr.distID][x.distID] + D[p.nbr.distID][y.distID]) / 4.0;
        	        }
        		}
        	}
        	return ((double) num_clusters - 1 - 2.0) * Dpu - clusterRowSum - (subtractedRowSum + Dpu);
        } else if (null == x.nbr) {     /* X is isolated,  Y  is not isolated*/
            return agg3wayLocal(x, y, y.nbr, testNode, subtractedRowSum, D, netNodes, num_clusters, num_active);
//            num_nodes += 2;
//            num_active--;
//            num_clusters--;
//            x.positionID = -1;
//            y.positionID = -1;
//            y.nbr.positionID = -1;
//            if (y.nbr != null) {y.nbr.positionID = -1;}
        } else if ((null == y.nbr) || (num_active == 4)) { /* Y is isolated,  X is not isolated
                                                    OR theres only four active nodes and none are isolated */
            return agg3wayLocal(y, x, x.nbr, testNode, subtractedRowSum, D, netNodes, num_clusters, num_active);
//            num_nodes += 2;
//            num_active--;
//            num_clusters--;
//            x.positionID = -1;
//            y.positionID = -1;
//            x.nbr.positionID = -1;
        } else {  /* Both nodes are connected to others and there are more than 4 active nodes */
            return agg4wayLocal(x.nbr, x, y, y.nbr, testNode, subtractedRowSum, D, netNodes, num_clusters, num_active);
//            num_nodes += 4;
//            num_active -= 2;
//            num_clusters--;
        }
	}
	
	private double agg3wayLocal(NetNode x, NetNode y, NetNode z, NetNode testNode, double subtractedRowSum, double[][] D, NetNode[] netNodes, int num_clusters, int num_active) {
        double clusterRowSum = 0;
        double clustDistance = 0;
		for (int i = 0; i < num_active; i++) { //Is this correct?
        	NetNode p = netNodes[i]; 
        	double myDist = 0;
    		if ((p.nbr == null || p.nbr.id > p.id) && (p != x && p != y && p != z)) {
    			double Dup = (2.0 / 3.0) * D[x.distID][p.distID] + D[y.distID][p.distID] / 3.0;
                double Dvp = (2.0 / 3.0) * D[z.distID][p.distID] + D[y.distID][p.distID] / 3.0;
    	        if (p.nbr == null) {// && (u.nbr != null))
    	        	myDist = (Dup + Dvp) / 2.0;
    	        } else {
    	        	double Duq = (2.0 / 3.0) * D[x.distID][p.nbr.distID] + D[y.distID][p.nbr.distID] / 3.0;
                    double Dvq = (2.0 / 3.0) * D[z.distID][p.nbr.distID] + D[y.distID][p.nbr.distID] / 3.0;
    	        	myDist = (Dup + Dvp + Duq + Dvq) / 4.0;
    	        }
    		}
    		if (myDist != 0) {
    			clusterRowSum += myDist;
    		}
    		if (p == testNode) {
    			clustDistance = myDist;
    		}
        }
		return ((double) num_clusters - 1 - 2.0) * clustDistance - clusterRowSum - (subtractedRowSum + clustDistance);
			
	}
	
	private double agg4wayLocal(NetNode x2, NetNode x, NetNode y, NetNode y2, NetNode testNode, 
			double subtractedRowSum, double[][] D, NetNode[] netNodes, int num_clusters, int num_active) {
//        u = agg3way(x2, x, y, amalgs, D, netNodes, num_nodes, num_active); /* Replace x2,x,y by two nodes, equal to x2_prev.next and y_prev.next. */
//        num_nodes += 2;
//        v = agg3way(u, u.nbr, y2, amalgs, D, netNodes, num_nodes, num_active-1); /* z = y_prev . next */
//        num_nodes += 2;
		double[] DupVec = new double[num_active];
		double[] DvpVec = new double[num_active];
		for (int i = 0; i < num_active; i++) {  //Subtract off the old value.  Add the new value.
        	NetNode p = netNodes[i];  //What about the 4-way? //What if p = u or v?
        	if (p != x) {
        		DupVec[i] = (2.0 / 3.0) * D[x2.distID][p.distID] + D[x.distID][p.distID] / 3.0;
        		DvpVec[i] = (2.0 / 3.0) * D[y.distID][p.distID] + D[x.distID][p.distID] / 3.0;
        	} else {
        		DupVec[i] = DvpVec[i] = 0;
        	}
        }
		//u1 = x2, v1 = y
		double[] DupVec2 = new double[num_active];
		double[] DvpVec2 = new double[num_active];
		double clusterRowSum = 0;
        double clustDistance = 0;
		for (int i = 0; i < num_active; i++) {
			NetNode p = netNodes[i];
			double myDist = 0;
			//(2.0 / 3.0) * D[u.distID][p.distID] + D[u.nbr.distID][p.distID] / 3.0;
			//(2.0 / 3.0) * D[y2.distID][p.distID] + D[u.nbr.distID][p.distID] / 3.0;
			if (p != x && p != y && p != y2 && p != x2) {
				if ((p.nbr == null || p.nbr.id > p.id)) {
					DupVec2[i] = (2.0 / 3.0) * DupVec[i] + DvpVec[i] / 3.0;
					DvpVec2[i] = (2.0 / 3.0) * D[y2.distID][p.distID] + DvpVec[i] / 3.0;
					if (p.nbr == null) {
						myDist = (DupVec2[i] + DvpVec2[i]) / 2.0;
					} else {
						int j = p.nbr.positionID;
						DupVec2[j] = (2.0 / 3.0) * DupVec[j] + DvpVec[j] / 3.0;
						DvpVec2[j] = (2.0 / 3.0) * D[y2.distID][p.nbr.distID] + DvpVec[j] / 3.0;
						myDist = (DupVec2[i] + DvpVec2[i] + DupVec2[j] + DvpVec2[j]) / 4.0; 
					}
				}
			} else {
				DupVec2[i] = DvpVec2[i] = 0;
			}
			clusterRowSum += myDist;
			if (p == testNode) {
				clustDistance = myDist;
			}
		}
		return ((double) num_clusters - 1 - 2.0) * clustDistance - clusterRowSum - (subtractedRowSum + clustDistance);
		
	}
	
	private double calculateClusterDistLocal(NetNode p, NetNode u) {
		Double Dpu = 0.0;
        if (p.nbr == null) {// && (u.nbr != null))
            Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID]) / 2.0;
        } else {
            Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID] + D[p.nbr.distID][u.distID] + D[p.nbr.distID][u.nbr.distID]) / 4.0;
        }
        return Dpu;
	}
}

//public static void shuffle(List<?> list, Random rnd) {
//    int size = list.size();
//    if (size < SHUFFLE_THRESHOLD || list instanceof RandomAccess) {
//        for (int i=size; i>1; i--)
//            swap(list, i-1, rnd.nextInt(i));
//    } else {
//        Object arr[] = list.toArray();
//
//        // Shuffle array
//        for (int i=size; i>1; i--)
//            swap(arr, i-1, rnd.nextInt(i));
//
//        // Dump array back into list
//        ListIterator it = list.listIterator();
//        for (int i=0; i<arr.length; i++) {
//            it.next();
//            it.set(arr[i]);
//        }
//    }
//}
//
//
// private static void swap(Object[] x, int a, int b) {
//    Object t = x[a];
//    x[a] = x[b];
//    x[b] = t;
//}