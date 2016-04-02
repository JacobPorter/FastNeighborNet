package nnet;

import java.util.Comparator;
import java.util.Stack;
import java.util.PriorityQueue;
import java.util.concurrent.ExecutorService;
import java.util.Arrays;

class NeighborNetFilter extends NetMakerOriginal {
	
	class BinTriple {
		final NetNode c1;
		final NetNode c2;
		final NetNode n1;
		final NetNode n2;
		final double distance;
		
		BinTriple(NetNode c1, NetNode c2, double distance) {
			this.c1 = c1;
			this.n1 = c1.nbr;
			this.c2 = c2;
			this.n2 = c2.nbr;
			this.distance = distance;
		}
		
		BinTriple() {
			this.c1 = null;
			this.c2 = null;
			this.n1 = null;
			this.n2 = null;
			this.distance = 0.0;
		}
		
		final boolean save() {
			return (c1.positionID >= 0 && c2.positionID >= 0 && c1.nbr == n1 && c2.nbr == n2 && c1.nbr != c2); // Is last one necessary?
		}
		
	}
	
	class BinHeap {
		final PriorityQueue<BinTriple> myQueue  = new PriorityQueue<BinTriple>(64, btc);
		double R1Max;
		double R2Max;
	}
	
	class BinHeaps {
		final int binAmount;
		final BinHeap heaps[][];
		double Rmax;
		double Rmin;
		double intervalLength;
		
		BinHeaps(final int binAmount) {
			this.binAmount = binAmount;
			heaps = new BinHeap[this.binAmount][];
			for (int i = 0; i < binAmount; i++) {
				heaps[i] = new BinHeap[i+1];
				for (int j = 0; j <= i; j++) {  //Correct??
					heaps[i][j] = new BinHeap();
				}
			}
		}
		
		void initialize(NetNode[] netNodes, int num_active) {
			Rmax = Double.MIN_VALUE;
			Rmin = Double.MAX_VALUE;
			//double i_rowD;
			//int clusters[] = new int[num_active];
//			double rowSums[] = new double[numTaxa]; //Pointless?
			//int clusterPointer = 0;
			for (int i = 0; i < num_active; i++) {
				NetNode p = netNodes[i];
				//clusterPointer++;
					if (p.Sx > Rmax) {
						Rmax = p.Sx;
					} else if (p.Sx < Rmin) {
						Rmin = p.Sx;
					}
			}
			intervalLength = (Rmax - Rmin) / binAmount;
			double R1current = Rmin;
			double R2current = Rmin;
			//System.err.println("Rmin and Rmax: " + Rmin + " " + Rmax);
			//can be multithreaded
			//Initialize and reset heaps
			for (int i = 0; i < binAmount; i++) {
				if (i!=0) {R1current += intervalLength;}
				R2current = Rmin + intervalLength * i;
				for (int j = 0; j  <= i; j++) {
					if (j!=i) {R2current += intervalLength;}
					heaps[i][j].myQueue.clear();
					heaps[i][j].R1Max = R1current + intervalLength;
					//heaps[i][j].R1Min = R1current;
					heaps[i][j].R2Max = R2current + intervalLength;
					//heaps[i][j].R2Min = R2current;
					// Is this necessary??
//					if (i == binAmount - 1) {
//						heaps[i][j].R1Max += eps;
//					}
//					if (j == binAmount - 1) {
//						heaps[i][j].R2Max += eps;
//					}
				}
			}
		}
		
		void insertBinTriple(BinTriple[] btarray) {
			for(BinTriple myBT: btarray) {
				if (myBT == null) {
					break;
				}
				double rowSum1 = myBT.c1.Sx;
				int index1 = (rowSum1 == Rmax) ? binAmount - 1 : (int) ((rowSum1 - Rmin) / intervalLength);
				index1 = index1 < 0 ? 0 : index1;
				double rowSum2 = myBT.c2.Sx;
				int index2 = (rowSum2 == Rmax) ? binAmount - 1 : (int) ((rowSum2 - Rmin) / intervalLength);
				index2 = index2 < 0 ? 0 : index2;
				//System.err.println(rowSum1 + " " + index1 + " " + rowSum2 + " " + index2);
				BinHeap myBinHeap = (index2>index1) ? heaps[index2][index1] : heaps[index1][index2]; 
				myBinHeap.myQueue.add(myBT);
			}
		}
	}
	
	class BinHeapsFreq extends BinHeaps {
		
		final double[] heapWidths;
		
		BinHeapsFreq(final int binAmount) {
			super(binAmount);
			heapWidths = new double[this.binAmount+1];
			//this.binAmount = binAmount;
		}
		
		@Override
		void initialize(NetNode[] netNodes, int num_active) {
			//Rmax = Double.MIN_VALUE;
			//Rmin = Double.MAX_VALUE;
			//Arrays.fill(heapWidths, null);
			Arrays.sort(netNodes, 0, num_active, null);
			Rmin = netNodes[0].Sx;
			Rmax = netNodes[num_active-1].Sx;
			int frequency = (int) num_active / binAmount;
			heapWidths[0] = Rmin;
			int counter = 1;
			for(int i = frequency; i < num_active && counter < heapWidths.length; i += frequency, counter++) {
				heapWidths[counter] = (netNodes[i-1].Sx + netNodes[i].Sx) / 2;
			}
			heapWidths[binAmount] = Rmax;
			System.out.println(Arrays.toString(heapWidths));
			System.out.println(frequency);
			int bacon = 0;
			for (NetNode nn : netNodes) {
				if (bacon % frequency == 0) {
					System.out.print(" |" + (int) bacon / frequency + "| ");
				}
				System.out.print(nn.Sx + " ");
				bacon++;
			}
			System.out.print("\n");
			for (int i = 0; i < binAmount; i++) {
				for (int j = 0; j  <= i; j++) {
					heaps[i][j].myQueue.clear();
					heaps[i][j].R1Max = heapWidths[i+1];
					heaps[i][j].R2Max = heapWidths[j+1];
				}	
			}
		}
		
		@Override
		void insertBinTriple(BinTriple[] btarray) {
			for(BinTriple myBT: btarray) {
				if (myBT == null) {
					break;
				}
				double rowSum1 = myBT.c1.Sx;
				double rowSum2 = myBT.c2.Sx;
				int index1 = findBucket(rowSum1);
				int index2 = findBucket(rowSum2);
				BinHeap myBinHeap = (index2>index1) ? heaps[index2][index1] : heaps[index1][index2]; 
				myBinHeap.myQueue.add(myBT);
			}
		}
		
		private int findBucket(double rowSum) {
			if (rowSum == Rmax) {
				return binAmount - 1;
			} else if (rowSum < Rmin) {
				return 0;
			}
			int low = 0;
			int middle = (int) binAmount / 2;
			int high = binAmount + 1;
			while (true) {
				if (low+1 == high || (rowSum >= heapWidths[low] && rowSum < heapWidths[low+1])) {
					return low;
				}
				double lowSum = heapWidths[low];
				double middleSum = heapWidths[middle];
				//double highSum = heapWidths[high];
				if (lowSum <= rowSum && rowSum < middleSum) {
					high = middle;
				} else {
					low = middle;
				}
				middle = (int) (low + high) / 2;
			}
		}
	}
	
	
	class BinTripleComparator implements Comparator<BinTriple> {
		public int compare(BinTriple left, BinTriple right) {
			if (left.distance < right.distance) {
				return -1;
			} else if (left.distance > right.distance) {
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	
	final BinTripleComparator btc = new BinTripleComparator();
	//final BinHeap heaps[][];
	final BinHeaps myBinHeaps;
	//Candidate list of BinTriples pulled off of the bin heaps
	BinTriple candidates[];
	// Size of candidate array
	int candidateSize = 1000;
	//Queue of positions in the candidates array that are available to put new BinTriples into
//	final ConcurrentLinkedQueue<Integer> availableCandidatePositions = new ConcurrentLinkedQueue<Integer>();
	//Farthest position in the candidates array to put a BinTriple into
	int nextCandidate = 0;
	//The number of iterations until the next rebuild
	int stepsUntilRebuild = 0;
	//Used to compute the ratio of iterations remaining for rebuilding the bin heaps
	double rebuildRatio = 0.5;
	//Track bin triple minimizer object
	BinTriple myMinBT;
	//Dummy bin triple object
	final BinTriple dummy = new BinTriple();
	//Block size for bin-heaps
	//final int binHeapBlockSize;
	//Class for coordinating work in finding valid clusters
//	final ValidClusters myVC;
	//double Rmin;
	//double Rmax;
	//double intervalLength;
	//Number of bins
	final int binAmount = 30;
	// Bin heap adder classes for parallelization
	//final addNewBTRunner[] myBTRunners;
	//Array representing iteration count of last update for clusterID
	final int[] binAssignments; //numTaxa
	//Array for multithreaded resetting of BinHeaps
	//final ResetBinHeaps[] myResetBinHeaps;
	//Array for adding clusters to the BinHeap
	//final AddClustersToBinHeap[] myAddClustersToBinHeap;
	//Array for finding the maximum and minimum row sum
	//final FindRowMaxMin[] myFindRowMaxMin;
	int iterCount;
	
	public NeighborNetFilter(double[][] d, int numTaxa, int numThreads, ExecutorService pool, boolean type) {
		super(d, numTaxa, numThreads, pool);
		//myVC = new ValidClusters(numTaxa);
		binAssignments = new int[numTaxa];
		//candidates = new BinTriple[candidateSize];
		if (type) {
			myBinHeaps = new BinHeaps(binAmount);
		} else {
			myBinHeaps = new BinHeapsFreq(binAmount);
		}
	}
		
	
	@Override
	protected void initialize(double D[][], NetNode[] netNodes, int num_nodes, int num_active, int num_clusters) {
		final BinTriple[] btlist = new BinTriple[(num_active * (num_active-1)) / 2];
		int index = 0;
    	for (int i = 0; i < num_active; i++) { //(NetNode p = netNodes.next; p != null; p = p.next) {
            NetNode p = netNodes[i]; 
        	if (p.nbr == null || p.nbr.id > p.id) {
        		for (int j = p.positionID+1; j < num_active; j++) {
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
                        if (netNodes.length == num_active) {// Check if this is the first time through or not.
                        	p.Sx += Dpq;
	                        if (p.nbr != null)
	                            p.nbr.Sx += Dpq;
	                        q.Sx += Dpq;
	                        if (q.nbr != null)
	                            q.nbr.Sx += Dpq;
                        }
                        //myBinHeaps.insertOneBinTriple(new BinTriple(p, q, Dpq));
                        btlist[index] = new BinTriple(p, q, Dpq);
                        index++;
                    }
                }
            }
        }
    	myBinHeaps.initialize(netNodes, num_active);
    	myBinHeaps.insertBinTriple(btlist);
    	stepsUntilRebuild = (int) (num_clusters*rebuildRatio);
    	candidates = new BinTriple[candidateSize];
    	nextCandidate = 0;
	}

	@Override
	protected void updateClusterDistances(NetNode u, double D[][], NetNode[] netNodes, int num_active) {
    	u.Sx = 0;
        u.nbr.Sx = 0;
        BinTriple[] btlist = new BinTriple[num_active];
        int index = 0;
        for (int i = 0; i < num_active; i++) {
        	NetNode p = netNodes[i];
        	if  ((p.nbr == null || p.nbr.id > p.id) && (u.nbr != p) && (u != p)) {
        		Double Dpu = 0.0;
                if (p.nbr == null) {// && (u.nbr != null))
                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID]) / 2.0;
                } else {
                    Dpu = (D[p.distID][u.distID] + D[p.distID][u.nbr.distID] + D[p.nbr.distID][u.distID] + D[p.nbr.distID][u.nbr.distID]) / 4.0;
                }
                btlist[index] = new BinTriple(u, p, Dpu);
                index++;
                p.Sx += Dpu;
                if (p.nbr != null) {p.nbr.Sx += Dpu; }  // Could remove addition and set to p.Sx
                u.Sx += Dpu;
        	}
        }
        u.nbr.Sx = u.Sx;
        myBinHeaps.insertBinTriple(btlist);
    }

	/**
	 * Overwrites Cx and Cy with the pair to combine.
	 */
	@Override
	protected void findNodes(Stack<NetNode> amalgs, double[][] D,
			NetNode[] netNodes, int num_nodes, int num_active, int num_clusters) {
		//System.err.println("Inside filter findNodes");
		best = Double.MAX_VALUE;
		Cx = Cy = null;
		// Check if heaps need to be recreated
		if (stepsUntilRebuild <= 0) {
			//System.err.println("Recreating heaps...");
			initialize(D, netNodes, num_nodes, num_active, num_clusters);
		}
		stepsUntilRebuild--;
		//Search through candidate list
		//boolean fromCand = true;
		double constant = ((double) num_clusters - 2.0);
		for (int i = 0; i < nextCandidate; i++) {
			BinTriple candidateTriple = candidates[i];
			if (candidateTriple.save()) {
                double Qpq = constant * candidateTriple.distance - candidateTriple.c1.Sx - candidateTriple.c2.Sx;
                //count++;
                /* Check if this is the best so far */
                if (Cx == null || Qpq < best) { // && (candidateTriple.c1.nbr != candidateTriple.c2)) {  // Last clause necessary
                    Cx = candidateTriple.c1;
                    Cy = candidateTriple.c2;
                    best = Qpq;
                }
			} else {
				candidates[i] = candidates[nextCandidate - 1];
				candidates[nextCandidate - 1] = null;
				nextCandidate--;
				i--;
			}
		}
		//System.err.println("");
		//Search through the bin heaps
		for (int i = 0; i < myBinHeaps.binAmount; i++) {
			for (int j = 0; j <= i; j++) {
				PriorityQueue<BinTriple> myQueue = myBinHeaps.heaps[i][j].myQueue;
				double maxRowSum = myBinHeaps.heaps[i][j].R1Max + myBinHeaps.heaps[i][j].R2Max;
				BinTriple bt = myQueue.poll();  //Should poll and then add it back on if needed?
				while (bt != null) {
					if (bt.save()) {
						double q = constant * bt.distance;
						double qLimit = q - maxRowSum;
		                if (Cx == null || (qLimit < best)) { // && (bt.c1.nbr != bt.c2)) {  // Last clause necessary //Need bin heap Rmax instead optimality
		                	double Qpq = q - bt.c1.Sx - bt.c2.Sx;
		                	if (Cx == null || (Qpq < best)) {
			                	//fromCand = false;
			                	Cx = bt.c1;
			                    Cy = bt.c2;
			                    best = Qpq;
			        		}
		                	//myQueue.poll();
		                    candidates[nextCandidate] = bt;
		                    nextCandidate++;
			                if (nextCandidate >= candidates.length) {
		        				BinTriple temp[] = new BinTriple[candidates.length * 2];
		        				System.arraycopy(candidates, 0, temp, 0, candidates.length);
		        				candidates = temp;
			                }
		                } else {
		                	myQueue.add(bt);
							break;
						}
					} else {
						//myQueue.poll();
					}
					bt = myQueue.poll();
				}
			}
		}
			//iterCount++;
			//System.err.println("From candidate array: " + fromCand);
	}
	
	
	
}