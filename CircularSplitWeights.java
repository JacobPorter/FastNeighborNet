package nnet;
//import com.sun.istack.internal.Nullable;

//import splits.algorithms.util.optimization.AbstractFunction;
//import splits.core.TaxaSet;
//import splits.nexus.Distances;
//import splits.nexus.Splits;



import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import nnet.DistancesAndNames;
/**
 * Given a circular ordering and a distance matrix,
 * computes the unconstrained or constrained least square weighted splits
 * <p/>
 * For all vectors, the canonical ordering of pairs is (0,1),(0,2),...,(0,n-1),(1,2),(1,3),...,(1,n-1), ...,(n-1,n)
 * <p/>
 * (i,j) ->  (2n - i -3)i/2 + j-1  .
 * <p/>
 * Increase i -> increase index by n-i-2
 * Decrease i -> decrease index by n-i-1.
 * <p/>
 * <p/>
 * x[i][j] is the split {i+1,i+2,...,j} | -------
 */

//TODO: Add parallelization

class NumericalException extends Exception {
	  public NumericalException() { super(); }
	  public NumericalException(String message) { super(message); }
	  public NumericalException(String message, Throwable cause) { super(message, cause); }
	  public NumericalException(Throwable cause) { super(cause); }
}

class SplitAndWeight {
	BitSet split;
	double weight;
}

public class CircularSplitWeights {
    /* Epsilon constant for the conjugate gradient algorithm */
    static final double CG_EPSILON = 1e-8;  // JSP Original
	//static final double CG_EPSILON = 1.0E-4D;
     //ToDo: Change way that variance options are passed.

    /**
     * Computes optimal least squares weights for a circular set of splits
     * @param ordering   Circular ordering of taxa
     * @param dist        Distance matrix
     * @param var          Method for variance
     * @param constrained    True if split weights are constrained to be positive.
     * @param cutoff         Do not include any splits with less than this weight.
     * @return              Set of splits weighted by least squares weights
     */
     public List<SplitAndWeight> getWeightedSplits(int[] ordering,
    		 DistancesAndNames dist, String var, boolean constrained, double cutoff, long maxIterations, double grad_cutoff, boolean useMax) {
        int ntax = dist.getNtax();
        //int npairs = (ntax * (ntax - 1)) / 2;

        //TODO: handle small number of taxa
        //Handle n=1,2 separately.
//        if (ntax == 1)
//            return new Splits(1);
//        if (ntax == 2) {
//            Splits smallSplits = new Splits(2);
//            float d_ij = (float) dist.get(ordering[1], ordering[2]);
//            if (d_ij > 0.0) {
//                TaxaSet a = new TaxaSet();
//                a.set(ordering[1]);
//                smallSplits.add(a, d_ij);
//            }
//            return smallSplits;
//        }

        /* Re-order taxa so that the ordering is 0,1,2,...,n-1 */
        //TODO: replace setupD with a distance matrix fetch so that half the memory will be used
        double[] d = setupD(dist, ordering);
        double[] v = setupV(dist, var, ordering);
        
        double[] x = getWeights(ntax, ordering, d, v, constrained, maxIterations, grad_cutoff, useMax);
        
        List<SplitAndWeight> splits = new ArrayList<SplitAndWeight>(ntax);
        int index = 0;
        for (int i = 0; i < ntax; i++) {
        	BitSet split = new BitSet(ntax);
        	for (int j = i + 1; j < ntax; j++) {
        		split.set(ordering[j]); //Should this be j+1?
        		if (x[index] > cutoff) {
        			SplitAndWeight saw = new SplitAndWeight();
        			saw.split = (BitSet) split.clone();
        			saw.weight = x[index];
        			splits.add(saw);
        		}
        		index++;
        	}
        }
        System.err.println("Number of splits: " + splits.size());
//        System.err.println(Arrays.toString(x));
//        int mycount = 0;
//        for (double ii: x) {
//        	if (ii > 0.0) {
//        		mycount++;
//        	}
//        }
//        System.err.println("Number bigger than zero: " + mycount);
        
        //return x;

        /* Construct the splits with the appropriate weights */
//        System.out.println("HELLO!!!");
//        Splits s_splits = new Splits(ntax);
//        index = 0;
//        for (int i = 0; i < ntax; i++) {
//            TaxaSet t = new TaxaSet();
//            for (int j = i + 1; j < ntax; j++) {
//                t.set(ordering[j]);
//                if (x[index] > cutoff)
//                    s_splits.add(t, (float) (x[index]));
//                index++;
//
//            }
//        }
        
//        Splits s_splits = new Splits(ntax);
//        index = 0;
//        for (int i = 0; i < ntax; i++) {
//            TaxaSet t = new TaxaSet();
//            for (int j = i + 1; j < ntax; j++) {
//                t.set(ordering[j + 1]);
//                if (x[index] > cutoff)
//                    s_splits.add(t, (float) (x[index]));
//                index++;
//
//            }
//        }
//        for (int i = 0; i < s_splits.getNsplits(); i++) {
//        	System.out.println(s_splits.get(i+1));
//        	System.out.println(s_splits.getWeight(i+1));
//        }
//        System.err.println("!!!" + s_splits.toString());
        return splits;
//        return splits;
        
    }
     
//     static double[] getWeights(int ntax, int[] ordering, double[] d, double[] v, boolean constrained) {
//    	 return getWeights(ntax, ordering, d, v, constrained, 0, -0.0000001, false);
//     }
     
     static double[] getWeights(int ntax, int[] ordering, double[] d, double[] v, boolean constrained, long maxIterations, double grad_cutoff, boolean useMax) {
    	 int npairs = (ntax * (ntax - 1)) / 2;
    	 double[] x = new double[npairs];

         if (!constrained || useMax && maxIterations == 1)
             CircularSplitWeights.runUnconstrainedLS(ntax, d, x);
         else // do constrained optimization
         {
             /* Initialise the weight matrix */
             double[] W = new double[npairs];
             for (int k = 0; k < npairs; k++) {
                 if (v[k] == 0.0)
                     W[k] = 10E10;
                 else
                     W[k] = 1.0 / v[k];
             }

             
                 /* Find the constrained optimal values for x */
                 runActiveConjugate(ntax, d, W, x, maxIterations, grad_cutoff, useMax);
//                 System.out.print("optimalx = [");
//                 for(int i=0;i<x.length;i++)
//                     System.out.print(" "+x[i]);
//                 System.out.println("]';");
         }
         return x;
    	 
     }
     
     

    /**
     * setup working distance so that ordering is trivial.
     * Note the the code assumes that taxa are labeled 0..ntax-1 and
     * we do the transition here. It is undone when extracting the splits
     *
     * @param dist     Distances block
     * @param ordering circular ordering
     * @return double[] distances stored as a vector
     */
    static double[] setupD(DistancesAndNames dist, int[] ordering) {
        return dist.get();
//    	int ntax = dist.getNtax();
//        double[] d = new double[(ntax * (ntax - 1)) / 2];
//        int index = 0;
//        for (int i = 0; i < ntax; i++)
//            for (int j = i + 1; j < ntax; j++)
//                d[index++] = dist.get(ordering[i], ordering[j]);
//        return d;
    }

    static double[] setupV(DistancesAndNames dist, String var, int[] ordering) {
        int ntax = dist.getNtax();
        int npairs = ((ntax - 1) * ntax) / 2;
        double[] v = new double[npairs];
        if (var.equalsIgnoreCase("ols")) {
        	Arrays.fill(v, 1.0);
        }

//        int index = 0;
//        for (int i = 0; i < ntax; i++)
//            for (int j = i + 1; j < ntax; j++) {
//                double dij = dist.get(ordering[i], ordering[j]);
//                //if (var.equalsIgnoreCase("ols"))
//                //    v[index] = 1.0;
//                else if (var.equalsIgnoreCase("fm1"))
//                    v[index] = dij;
//                else if (var.equalsIgnoreCase("fm2"))
//                    v[index] = dij * dij;
//                else
//                    v[index] = dist.getVar(ordering[i + 1], ordering[j + 1]);
//                index++;
//            }
        return v;
    }


    /**
     * Compute the branch lengths for unconstrained least squares using
     * the formula of Chepoi and Fichet (this takes O(N^2) time only!).
     *
     * @param n the number of taxa
     * @param d the distance matrix
     * @param x the split weights
     */
    static void runUnconstrainedLS(int n, double[] d, double[] x) {
        int index = 0;

        for (int i = 0; i <= n - 3; i++) {
            //index = (i,i+1)
            //x[i,i+1] = (d[i][i+1] + d[i+1][i+2] - d[i,i+2])/2
            x[index] = (d[index] + d[index + (n - i - 2) + 1] - d[index + 1]) / 2.0;
            index++;
            for (int j = i + 2; j <= n - 2; j++) {
                //x[i][j] = ( d[i,j] + d[i+1,j+1] - d[i,j+1] - d[i+1][j])
                x[index] = (d[index] + d[index + (n - i - 2) + 1] - d[index + 1] - d[index + (n - i - 2)]) / 2.0;
                index++;
            }
            //index = (i,n-1)

            if (i == 0) //(0,n-1)
                x[index] = (d[0] + d[n - 2] - d[2 * n - 4]) / 2.0; //(d[0,1] + d[0,n-1] - d[1,n-1])/2
            else
                //x[i][n-1] == (d[i,n-1] + d[i+1,0] - d[i,0] - d[i+1,n-1])
                x[index] = (d[index] + d[i] - d[i - 1] - d[index + (n - i - 2)]) / 2.0;
            index++;
        }
        //index = (n-2,n-1)
        x[index] = (d[index] + d[n - 2] - d[n - 3]) / 2.0;
    }

    /**
     * Returns the array indices for the smallest propKept proportion of negative values in x.
     * In the case of ties, priority is given to the earliest entries.
     * Size of resulting array will be propKept * (number of negative entries) rounded up.
     *
     * @param x        returns an array
     * @param propKept the
     * @return int[] array of indices
     */
    static int[] worstIndices(double[] x, double propKept) {


        if (propKept == 0)
            return null;

        int n = x.length;

        int numNeg = 0;
        for (int i = 0; i < n; i++)
            if (x[i] < 0.0) 
                numNeg++;

        if (numNeg == 0)
            return null;

        //Make a copy of negative values in x.
        double[] xcopy = new double[numNeg];
        int j = 0;
        for (int i = 0; i < n; i++)
            if (x[i] < 0.0)  
                xcopy[j++] = x[i];

        //Sort the copy
        Arrays.sort(xcopy);
        
        //Find the cut-off value. All values greater than this should
        //be returned, as well as some (or perhaps all) of the values
        //equal to this.
        int nkept = (int) Math.ceil(propKept * numNeg);  //Ranges from 1 to n
        double cutoff = xcopy[nkept - 1];
        
        //we now fill the result vector. Values < cutoff are filled
        //in from the front. Values == cutoff are filled in the back.
        //Values filled in from the back can be overwritten by values
        //filled in from the front, but not vice versa.
        int[] result = new int[nkept];
        int front = 0, back = nkept - 1;

        for (int i = 0; i < n; i++) {
            if (x[i] < cutoff)
                result[front++] = i; //Definitely in the top entries.
            else if (x[i] == cutoff) {
                if (back >= front)
                    result[back--] = i;
            }
        }
        return result;
    }


    static void printvec(String msg, double[] x) {
        int n = x.length;
        DecimalFormat fmt = new DecimalFormat("#0.00000");

        System.out.print(msg + "\t");
        for (int i = 0; i < n; i++)
            System.out.print(" " + fmt.format(x[i]));
        System.out.println();
    }


    /**
     * Uses an active set method with the conjugate gradient algorithm to find x that minimises
     * <p/>
     * (Ax - d)W(Ax-d)
     * <p/>
     * Here, A is the design matrix for the set of cyclic splits with ordering 0,1,2,...,n-1
     * d is the distance vector, with pairs in order (0,1),(0,2),...,(0,n-1),(1,2),(1,3),...,(1,n-1), ...,(n-1,n)
     * W is a vector of variances for d, with pairs in same order as d.
     * x is a vector of split weights, with pairs in same order as d. The split (i,j), for i<j, is {i,i+1,...,j-1}| rest
     *
     * @param ntax The number of taxa
     * @param d    the distance matrix
     * @param W    the weight matrix
     * @param x    the split weights
     */
    static void runActiveConjugate(int ntax, double[] d, double[] W, double[] x, long maxIterations, double grad_cutoff, boolean useMax) {
        final boolean collapse_many_negs = true;

        int npairs = d.length;
        if (W.length != npairs || x.length != npairs)
            throw new IllegalArgumentException("Vectors d,W,x have different dimensions");

        /* First evaluate the unconstrained optima. If this is feasible then we don't have to do anything more! */
        CircularSplitWeights.runUnconstrainedLS(ntax, d, x);
        boolean all_positive = true;
        for (int k = 0; k < npairs && all_positive; k++)
            if (x[k] < 0.0)
                all_positive = false;

        if (all_positive) /* If the unconstrained optimum is feasible then it is also the constrained optimum */
            return;

        /* Allocate memory for the "utility" vectors */
        double[] r = new double[npairs];
        double[] w = new double[npairs];
        double[] p = new double[npairs];
        double[] y = new double[npairs];
        double[] old_x = new double[npairs];
        Arrays.fill(old_x, 1.0);

        /* Initialise active - originally no variables are active (held to 0.0) */
        boolean[] active = new boolean[npairs];
        Arrays.fill(active, false);

        /* Allocate and compute AtWd */
        double[] AtWd = new double[npairs];
        for (int k = 0; k < npairs; k++)
            y[k] = W[k] * d[k];
        CircularSplitWeights.calculateAtx(ntax, y, AtWd);
        
        boolean first_pass = true; //This is the first time through the loops.
        maxIterations--;
        while (true) {
        	if (useMax) {
        		//System.err.print(maxIterations + " ");
        		maxIterations--;
        	}
        	//JSP: Use with debug switch
            if (false) {
            	System.err.println("Outer: nactive = " + countNumActive(active)+ " of "+npairs);
            }
            
            while (true) /* Inner loop: find the next feasible optimum */ {

            	//JSP: Use with debug code
            	if (false) {
            		System.err.println("\t Inner: nactive = " + countNumActive(active)+ " of "+npairs);
            	}


                if (!first_pass)  /* The first time through we use the unconstrained branch lengths */
                    CircularSplitWeights.circularConjugateGrads(ntax, npairs, r, w, p, y, W, AtWd, active, x);
                first_pass = false;
               
                if (collapse_many_negs) { /* Typically, a large number of edges are negative, so on the first
                								pass of the algorithm we add the worst 60% to the active set */
                    int[] entriesToContract = worstIndices(x, 0.6);
//                    for (int index = 0; index < npairs; index++) {  //JSP: added this loop to add indices with weight identically 0 to the active set 
//                    	if (x[index] == 0.0) {
//                    		active[index] = true;
//                    	}
//                    }
                    if (entriesToContract != null) {
                        int numToContract = entriesToContract.length;
                        for (int k = 0; k < numToContract; k++) {
                            int index =  entriesToContract[k];
                            x[index] = 0.0;
                            active[index] = true;
                        }
                        CircularSplitWeights.circularConjugateGrads(ntax, npairs, r, w, p, y, W, AtWd, active, x); /* Re-optimise, so that the current x is always optimal */
                    }
                }
                
                
                int min_i = -1;
                double min_xi = -1.0;
                for (int i = 0; i < npairs; i++) {
                    if (x[i] < 0.0) {
                        double xi = (old_x[i]) / (old_x[i] - x[i]);
                        if ((min_i == -1) || (xi < min_xi)) {
                            min_i = i;
                            min_xi = xi;
                        }
                    }
                }

                if (min_i == -1) /* This is a feasible solution - go to the next stage to check if its also optimal */
                    break;
                else {/* There are still negative edges. We move to the feasible point that is closest to
                							x on the line from x to old_x */

                    for (int i = 0; i < npairs; i++) /* Move to the last feasible solution on the path from old_x to x */
                        if (!active[i])
                            old_x[i] += min_xi * (x[i] - old_x[i]);
                    active[min_i] = true; /* Add the first constraint met to the active set */
                    x[min_i] = 0.0; /* This fixes problems with round-off errors */
                }
            }  /* End inner while loop */

            /* Find i,j that minimizes the gradient over all i,j in the active set. Note that grad = 2(AtWAb-AtWd) */
            calculateAb(ntax, x, y);
            for (int i = 0; i < npairs; i++)
                y[i] *= W[i];
            calculateAtx(ntax, y, r); /* r = AtWAx */

            /* We check to see that we are at a constrained minimum.... that is that the gradient is positive for
            * all i,j in the active set.
            */

            //double pgnorm = 0.0;

            int min_i = -1;
            double min_grad = 1.0;  // JSP should this be DOUBLE.MAX ?
            for (int i = 0; i < npairs; i++) {
                r[i] -= AtWd[i];
                r[i] *= 2.0;
                //if (!active[i] && x[i] < CG_EPSILON && x[i] > -CG_EPSILON ) {  //JSP: Added this condition
                //	active[i] = true;
                //}
                if (active[i]) {  // (active[i]) {  JSP: This is the original condition
                    double grad_ij = r[i];
                    if ((min_i == -1) || (grad_ij < min_grad)) {
                        min_i = i;
                        min_grad = grad_ij;
                    }
                }
                //if (!active[i] || r[i] < 0)
                //    pgnorm += r[i]*r[i];
            }
            
/*            boolean activeCheck = true;
            boolean inactiveCheck = true;
            for (int i = 0; i < npairs; i++) {
            	if (x[i] < CG_EPSILON && x[i] > -CG_EPSILON ) {
        			if (r[i] < -CG_EPSILON) {
        				activeCheck = false;
        			}
            	} else {
            		if (r[i] > CG_EPSILON || r[i] < -CG_EPSILON) {
        				inactiveCheck = false;
        			}
            	}
            }*/

            if ((min_i == -1) || (min_grad > -0.0000001)) { //(min_grad > -0.0000001)) { //|| (useMax && (maxIterations <= 0))) {  // ((min_i == -1) || (min_grad > -0.0000001))      {
            	System.err.println("min_i: " + min_i + " min_grad: " + min_grad);
            	double pgnorm = 0.0;
            	StringBuffer sb2000 = new StringBuffer("");
            	for (int i = 0; i < npairs; i++) {
            		if (x[i] < CG_EPSILON && x[i] > -CG_EPSILON ) {
            			if (r[i] < -CG_EPSILON) {
            				System.err.println("The gradient should be positive: " + i + " " + x[i] + " " + r[i]);
            			}
            		} else {
            			if (r[i] > CG_EPSILON || r[i] < -CG_EPSILON) {
            				System.err.println("The gradient should be zero: " + i + " " + x[i] + " " + r[i]);
            			}
            		}
            	}
            	for (int i = 0; i < npairs; i++) {
            		if (!active[i] || r[i] < 0)
                        pgnorm += r[i]*r[i];
            		if (active[i]) {
            			//System.err.print(i + ":" + r[i] + " ");
            			if (r[i] < CG_EPSILON ) {
            				sb2000.append(i + " ");
            			}
            		}
            	}
            	System.err.print("\n");
            	int activeCount = 0;
            	int inactiveCount = 0;
            	for (int i = 0; i < npairs; i++) {
            		if (active[i]) {
            			activeCount++;
            		} else {
            			inactiveCount++;
            		}
            		if (active[i] && (x[i] > CG_EPSILON || x[i] < -CG_EPSILON )) {
            			System.err.println("Non zero weight is active: " + i + " " + x[i]);
            		}
            		if (!active[i] && (x[i] < CG_EPSILON && x[i] > -CG_EPSILON ) ) {
            			System.err.println("Zero weight is not active: " + i + " " + x[i] + " " + r[i]);
            		}
            	}
            	System.err.println("Locations of negative gradients that should not be negative: " + sb2000);
            	System.err.println("pg norm is "+pgnorm);
            	System.err.println("Number of active: " + activeCount + " number of inactive: " + inactiveCount);
            	return; /* We have arrived at the constrained optimum */
            }
            else
                active[min_i] = false;
        }
    }

    static int countNumActive(boolean[] active) {
        int n=active.length;
        int c = 0;
        for (int i=0;i<n;i++)
            if (active[i])
                c++;
        return c;
    }


/* Compute the row sum in d. */

    static double rowsum(int n, double[] d, int k) {
        double r = 0;
        int index = 0;

        if (k > 0) {
            index = k - 1;     //The index for (0,k)

            //First sum the pairs (i,k) for i<k
            for (int i = 0; i < k; i++) {
                r += d[index];
                index += (n - i - 2);
            }
            index++;
        }
        //we now have index = (k,k+1)
        //Now sum the pairs (k,j) for k<j
        for (int j = k + 1; j < n; j++)
            r += d[index++];

        return r;
    }


    /**
     * Computes p = A^Td, where A is the topological matrix for the
     * splits with circular ordering 0,1,2,....,ntax-1
     * *
     *
     * @param n number of taxa
     * @param d distance matrix
     * @param p the result
     */
    static void calculateAtx(int n, double[] d, double[] p) {

//First the trivial splits
        int index = 0;
        for (int i = 0; i < n - 1; i++) {
            p[index] = rowsum(n, d, i + 1);
            index += (n - i - 1);
        }

        //Now the splits separating out two.
        index = 1;
        for (int i = 0; i < n - 2; i++) {
            //index = (i,i+2)

            //p[i][i+2] = p[i][i+1] + p[i + 1][i + 2] - 2 * d[i + 1][i + 2];
            p[index] = p[index - 1] + p[index + (n - i - 2)] - 2 * d[index + (n - i - 2)];
            index += (n - i - 2) + 1;
        }

        //Now the remaining splits
        for (int k = 3; k <= n - 1; k++) {
            index = k - 1;
            for (int i = 0; i <= n - k - 1; i++) {
                //index = (i,i+k)

                // p[i][j] = p[i][j - 1] + p[i+1][j] - p[i+1][j - 1] - 2.0 * d[i+1][j];
                p[index] = p[index - 1] + p[index + n - i - 2] - p[index + n - i - 3] - 2.0 * d[index + n - i - 2];
                index += (n - i - 2) + 1;
            }
        }
    }

    /**
     * Computes d = Ab, where A is the topological matrix for the
     * splits with circular ordering 0,1,2,....,ntax-1
     *
     * @param n number of taxa
     * @param b split weights
     * @param d pairwise distances from split weights
     */
    static void calculateAb(int n, double[] b, double[] d) {
        double d_ij;

        //First the pairs distance one apart.
        int index;
        int dindex = 0;

//        HashMap<Integer, Double> indexesDone = new HashMap<Integer, Double>();
        for (int i = 0; i <= n - 2; i++) {
            d_ij = 0.0;
            //Sum over splits (k,i) 0<=k<i.
            index = i - 1;  //(0,i)
            for (int k = 0; k <= i - 1; k++) {
                d_ij += b[index];  //(k,i)
                index += (n - k - 2);
            }
            index++;
            //index = (i,i+1)
            for (int k = i + 1; k <= n - 1; k++)  //sum over splits (i,k)  i+1<=k<=n-1
                d_ij += b[index++];

            d[dindex] = d_ij;
//            indexesDone.put(dindex, d_ij);
            dindex += (n - i - 2) + 1;
        }

        //Distances two apart.
        index = 1; //(0,2)
        for (int i = 0; i <= n - 3; i++) {
//            d[i ][i+2] = d[i ][i+1] + d[i + 1][i + 2] - 2 * b[i][i+1];

            d[index] = d[index - 1] + d[index + (n - i - 2)] - 2 * b[index - 1];
//            indexesDone.put(dindex, d[index]);
            index += 1 + (n - i - 2);
        }

        //return;
//        try {
//        FileWriter fw = new FileWriter(new File("TestCSW_Ab.txt"));
//		BufferedWriter bw = new BufferedWriter(fw);
//        bw.write("Testing Ab...\r\n");
//        int count = 0;
//        HashMap<Integer, Double> thisIteration = new HashMap<Integer, Double>();
//        String notSeen = "";
        for (int k = 3; k <= n - 1; k++) {
            index = k - 1;
            for (int i = 0; i <= n - k - 1; i++) {
                //int j = i + k;
                //d[i][j] = d[i][j - 1] + d[i+1][j] - d[i+1][j - 1] - 2.0 * b[i][j - 1];
//            	String print_indices = "";
//            	
//            	print_indices += " IND:" + (index);
//            	if (thisIteration.containsKey(index-1)) {
//            		//printStuff = true;
//            		print_indices += " " + (index-1);
//            		count++;
//            	} else {
//            		print_indices += " !" + (index-1);
//            	}
//            	if (thisIteration.containsKey(index + (n - i - 2))) {
//            		//printStuff = true;
//            		print_indices += " " + (index + (n - i - 2));
//            		count++;
//            	} else {
//            		print_indices += " !" + (index + (n - i - 2));
//            	}
//            	if (thisIteration.containsKey(index + (n - i - 2) - 1)) {
//            		//printStuff = true;
//            		print_indices += " " + (index + (n - i - 2) - 1);
//            		count++;
//            	} else {
//            		print_indices += " !" + (index + (n - i - 2) - 1);
//            	} 
//            	
//            	if (true) {
//            		bw.write(print_indices + "\r\n");
//            	}
                d[index] = d[index - 1] + d[index + (n - i - 2)] - d[index + (n - i - 2) - 1] - 2.0 * b[index - 1];
//                thisIteration.put(index, d[index]);
                index += 1 + (n - i - 2);
            }
        }
//        bw.write("\r\n");
//        bw.write(count + "\r\n");
//        bw.flush();
//        bw.close();
//        } catch (Exception e) {}
        
    }


    /**
     * Computes sum of squares of the lower triangle of the matrix x
     *
     * @param x the matrix
     * @return sum of squares of the lower triangle
     */
    static double norm(double[] x) {
        int n = x.length;
        double ss = 0.0;
        double xk;
        for (int k = 0; k < n; k++) {
            xk = x[k];
            ss += xk * xk;
        }
        return ss;
    }


    /**
     * Conjugate gradient algorithm solving A^tWA x = b (where b = AtWd)
     * such that all x[i][j] for which active[i][j] = true are set to zero.
     * We assume that x[i][j] is zero for all active i,j, and use the given
     * values for x as our starting vector.
     *
     * @param ntax   the number of taxa
     * @param npairs dimension of b and x
     * @param r      stratch matrix
     * @param w      stratch matrix
     * @param p      stratch matrix
     * @param y      stratch matrix
     * @param W      the W matrix
     * @param b      the b matrix
     * @param active the active constraints
     * @param x      the x matrix
     */
    static void circularConjugateGrads(int ntax, int npairs,
                                               double[] r, double[] w, double[] p, double[] y,
                                               double[] W, double[] b,
                                               boolean[] active, double[] x) {
        int kmax = ntax * (ntax - 1) / 2;
/* Maximum number of iterations of the cg algorithm (probably too many) */

        calculateAb(ntax, x, y);

        for (int k = 0; k < npairs; k++)
            y[k] = W[k] * y[k];
        calculateAtx(ntax, y, r); /*r = AtWAx */

        for (int k = 0; k < npairs; k++)
            if (!active[k])
                r[k] = b[k] - r[k];
            else
                r[k] = 0.0;

        double rho = norm(r);
        double rho_old = 0;

        double e_0 = CG_EPSILON * Math.sqrt(norm(b));
        int k = 0;

        while ((rho > e_0 * e_0) && (k < kmax)) {

            k = k + 1;
            if (k == 1) {
                for (int i = 0; i < npairs; i++)
                    p[i] = r[i];

            } else {
                double beta = rho / rho_old;
                //System.out.println("bbeta = " + beta);
                for (int i = 0; i < npairs; i++)
                    p[i] = r[i] + beta * p[i];

            }

            calculateAb(ntax, p, y);
            for (int i = 0; i < npairs; i++)
                y[i] *= W[i];

            calculateAtx(ntax, y, w); /*w = AtWAp */
            for (int i = 0; i < npairs; i++)
                if (active[i])
                    w[i] = 0.0;

            double alpha = 0.0;
            for (int i = 0; i < npairs; i++)
                alpha += p[i] * w[i];
            alpha = rho / alpha;

/* Update x and the residual, r */
            for (int i = 0; i < npairs; i++) {
                x[i] += alpha * p[i];
                r[i] -= alpha * w[i];
            }
            rho_old = rho;
            rho = norm(r);
        }
    }


    /**
     * Implements a function efficiently evaluating the sum of squares when fitting circular
     * least squares.
     */
//    public class CircularLSfunction { //extends AbstractFunction  {
//
//        private double[] d; //Re-ordered distance matrix
//        private double[] W; //Re-ordered weight matrix
//        private int ntax;
//        private int npairs;
//        private double[] scratch;
//        private double[] atwd; //AtWd
//        private double[] xs; //x shifted to start at 1.
//
//        public CircularLSfunction(int ntax, double[] d, double[] W) {
//            this.ntax = ntax;
//            this.npairs = ntax*(ntax-1)/2;
//            this.d = d;
//            this.W = W; //Should we be making copies of these?
//            this.scratch = new double[npairs];
//            this.atwd = new double[npairs];
//            for(int k=0;k<npairs;k++)
//                scratch[k] = W[k]*d[k];
//            calculateAtx(ntax,scratch,atwd);
//            xs = new double[npairs];
//        }
//
//        /**
//         * Computes xtAtWAx - 2xtAtWd, where A is the topological matrix for the circular splits.
//         * @param x  array of doubles
//         * @return    double
//         */
//        //@Override
//        public double get_val(double[] x) {
//
//            System.arraycopy(x,1,xs,0,npairs);
//            calculateAb(ntax,xs,scratch);
//            double xtAtWAx = 0.0;
//            double xtAtWd = 0.0;
//            for(int k=0;k<npairs;k++)
//                xtAtWAx += scratch[k] * scratch[k] * W[k];
//
//            for(int k=0;k<npairs;k++)
//                xtAtWd += xs[k]*atwd[k];
//
//            return xtAtWAx - 2.0*xtAtWd;
//        }
//
//        /**
//         * Computes
//         * 2AtWAx - 2AtWd
//         * @param x  1d array of doubles
//         * @param g  1d array with the same dimension as x
//         *
//         */
//        //@Override
//        public void get_grad(double[] x, double[] g) {
//            System.arraycopy(x,1,xs,0,npairs);
//
//            calculateAb(ntax,xs,scratch);
//            for(int k=0;k<npairs;k++)
//                scratch[k] *= W[k];
//            calculateAtx(ntax,scratch,g);
//            for(int k=0;k<npairs;k++)
//                g[k+1] = 2*(g[k] - atwd[k]);
//            g[0] = 0.0;
//        }
//
//
//        /**
//         * Same as the above, except we avoid one call to calculateAb
//         * @param x   1d array of double
//         * @param g   array of double with same dimension as x. Overwritten with the gradient
//         * @return    value of function at x
//         */
//        //@Override
//        public double get_val_and_grad(double[] x, double[] g) {
//            System.arraycopy(x,1,xs,0,npairs);
//
//
//            calculateAb(ntax,xs,scratch);
//            double xtAtWAx = 0.0;
//            double xtAtWd = 0.0;
//            for(int k=0;k<npairs;k++)
//                xtAtWAx += scratch[k] * scratch[k] * W[k];
//            for(int k=0;k<npairs;k++)
//                xtAtWd += xs[k]*atwd[k];
//
//            for(int k=0;k<npairs;k++)
//                scratch[k] *= W[k];
//            calculateAtx(ntax,scratch,g);
//            for(int k=0;k<npairs;k++)
//                g[k+1] = 2*(g[k] - atwd[k]);
//            g[0] = 0.0;
//
//            return xtAtWAx - 2.0*xtAtWd;
//        }
//
//        /**
//         * Returns 2A'WAv
//         * @param x 1d array of double
//         * @param v  1d array of double with the same dimensions as x
//         * @param hv  1d array of double with the same dimensions as x
//         *
//         * The array hv is overwritten with the produce H(x)v, where H(x) is the Hessian
//         */
//        //@Override
//        public void get_Hv(@Nullable double[] x, double[] v, double[] hv) {
//
//            double[] vs = xs;
//            System.arraycopy(v,1,vs,0,npairs);
//
//            calculateAb(ntax,vs,scratch);
//            for(int k=0;k<npairs;k++)
//                scratch[k] *= W[k];
//            calculateAtx(ntax,scratch,hv);
//
//            for(int k=0;k<npairs;k++)
//                hv[k+1]=hv[k]*2.0;
//            hv[0] = 0;
//        }
//
//        public void printMatlabDebug()  {
//            int n = ntax;
//            int npairs = ntax*(n-1)/2;
//            int index1=0,index2=0;
//
//
//            //Print the X matrix.
//            System.out.println("X = [");
//            for (int i=0;i<n;i++) {
//                for(int j=i+1;j<n;j++) {
//                    //Split (i+1,...,j).
//                    for (int k=0;k<n;k++) {
//                        for(int l=k+1;l<n;l++) {
//                            if ((k<=i && l>i && l<= j)||(i<k && k<=j && j<l))
//                                System.out.print("1 ");
//                            else
//                                System.out.print("0 ");
//                        }
//                    }
//                    if (i<n-2)
//                        System.out.println();
//                    else
//                        System.out.println("]';");
//                }
//            }
//
//            //Print the w vector
//            System.out.print("w = [");
//            for(int i=0;i<npairs;i++)
//                System.out.print(" "+W[i]);
//            System.out.println("]';");
//            System.out.println("W = diag(w);");
//
//            //Print the d vector
//             System.out.print("d = [");
//            for(int i=0;i<npairs;i++)
//                System.out.print(" "+d[i]);
//            System.out.println("]';");
//
//        }
//
//    }



}
