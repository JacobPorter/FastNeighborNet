package nnet;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;

//import cern.colt.Arrays;
//import java.io.BufferedWriter;

final public class DistancesAndNames {
	final double[] distances;
	final int nTaxa;
	final String[] names;
	//double[][] distances_2D;
	
	public DistancesAndNames(final double[] d, final String[] s, final int nTaxa) {
		distances = d;
		this.nTaxa = nTaxa;
		names = s;
	}
	
	public int upperIndex(final int i, final int j) {
		final int a;
		final int b;
		if (i == j) {
			return -1;
		}
		if (j < i) {
			a = j;
			b = i;
		} else {
			a = i;
			b = j;
		}
		return a*(nTaxa-1) - a*(a-1)/2 + b - (a+1);
	}
	
	
	
	
	public DistancesAndNames(final String fileString, final int numTaxa) throws Exception {
		this.nTaxa = numTaxa;
		final int npairs = (numTaxa * (numTaxa - 1)) / 2;
		distances = new double[npairs];
		//distances_2D = new double[nTaxa][];
		names = new String[nTaxa];
		//System.out.println("Number of taxa: " + numTaxa);
		//ComputeDistance[] myCDs = new ComputeDistance[numThreads];
		//List<Future> futures = null;
		try {
			Charset charset = Charset.forName("US-ASCII");
			try (BufferedReader reader = Files.newBufferedReader(Paths.get(fileString), charset)) {
			   String line = reader.readLine();
			   int row = 0;
			   //int count = 0;
			   String[] copy = new String[nTaxa];
			   //int currentBlockSize = 1, bi = 1, bj = 0, numBlock = 0;
			   //if (multiThread) {
			   // 	futures = new ArrayList<Future>(numThreads);
			   // }
			    while ((line = reader.readLine()) != null) {
			    	//double[] myRow = new double[row];
			    	String[] ss = line.split(" ");
			    	//System.out.println(Arrays.toString(ss));
			    	try {
			    		names[row] = ss[0];
			    	} catch (ArrayIndexOutOfBoundsException ai) {
			    		break;
			    	}
			    	int innerCount = 0;
			    	for (int i = 1 ; i < ss.length; i++ ) {
			    		if (!ss[i].trim().isEmpty()) {
			    			String[] temp = ss[i].split("\t");
			    			for (int j = 0 ; j < temp.length; j++) {
			    				copy[innerCount] = temp[j];
			    				innerCount++;
			    			}
			    		}
			    	}
			    	int column = 0;
			    	for ( ; column < row; column++) {
			    		//myRow[column] = Double.valueOf(copy[column]);
			    		int index = upperIndex(row, column);
			    		if (index > -1) {
			    			distances[index] = Double.valueOf(copy[column]);
			    		} else {
			    			System.err.println("Encountered a -1 for distances: " + row + ", " + column);
			    		}
//			    		if (multiThread) {
//			    			currentBlockSize++;
//			    			if (currentBlockSize == blockSize) {
//			    				myCDs[numBlock] = new ComputeDistance(bi, bj, row, column);
//			    				if (bi == bj) {
//			    					System.err.println("Error: block row, column the same: " + row);
//			    				}
//			    				bi = row; bj = column;
//			    				//System.err.println("numBlock " + numBlock + "bi: " + bi + " bj: " + bj);
//			    				currentBlockSize = 1; numBlock++;
//			    			}
//			    		}
			    		//count++;
			    	}
			    	//distances_2D[row] = myRow;
			    	row++;
			    }
			    //System.err.println("numBlock " + numBlock + "bi: " + row + " bj: " + numTaxa);
//			    if (numBlock != numThreads) {
//			    	myCDs[numBlock] = new ComputeDistance(bi, bj, row - 1, numTaxa-2);
//			    }
			} catch (IOException x) {
			    System.err.format("IOException: %s%n", x);
			}
		} finally {
		}
//		Charset charset = Charset.forName("US-ASCII");
//		try (BufferedWriter writer = Files.newBufferedWriter(Paths.get("Test_Distances.txt"), charset)) {
//			for (double D : distances) {
//				writer.write(D + " ");
//			}
//			writer.write("\n");
//			writer.flush();
//		}
//		for (int i = 0; i < nTaxa; i++) {
//			for (int j = i + 1; j < nTaxa; j++) {
//				if (distances[upperIndex(i, j)] != get(i, j)) {
//					System.err.println("Error getting distance at " + i + " " + j + " and " + upperIndex(i, j));
//				}
//			}
//		}
	}
	
	public int getNtax() {
		return this.nTaxa;
	}
	
	public double[] get() {
		return this.distances;
	}
	
	public double get(final int i, final int j) {
		final int uIndex = upperIndex(i, j);
		if (uIndex < 0) {
			return 0.0;
		} else {
			return distances[uIndex];
		}
	}
//	
//	public double getVar(int i, int j) {
//		return 0.0;
//	}
}