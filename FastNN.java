package nnet;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Locale;
import java.util.Arrays;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.ParseException;

import nnet.CircularSplitWeights;
import edu.rit.numeric.NonNegativeLeastSquares;

class NamesAndDistances {
	final double[][] distances;
	final String[] names;
	NamesAndDistances(String[] n, double[][] d) {
		distances = d;
		names = n;
	}
}

public class FastNN {
		
	protected static NamesAndDistances getFile(String fileString, int numTaxa) throws Exception {
		double[][] distances = new double[numTaxa][];
		String[] names = new String[numTaxa];
		//System.out.println("Number of taxa: " + numTaxa);
		//ComputeDistance[] myCDs = new ComputeDistance[numThreads];
		//List<Future> futures = null;
		try {
			Charset charset = Charset.forName("US-ASCII");
			try (BufferedReader reader = Files.newBufferedReader(Paths.get(fileString), charset)) {
			   String line = reader.readLine();
			   int row = 0;
			   //int count = 0;
			   String[] copy = new String[numTaxa];
			   //int currentBlockSize = 1, bi = 1, bj = 0, numBlock = 0;
			   //if (multiThread) {
			   // 	futures = new ArrayList<Future>(numThreads);
			   // }
			    while ((line = reader.readLine()) != null) {
			    	double[] myRow = new double[row];
			    	String[] ss = line.split(" ");
			    	names[row] = ss[0];
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
			    		myRow[column] = Double.valueOf(copy[column]);
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
			    	distances[row] = myRow;
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
		return new NamesAndDistances(names, distances);
	}
	

	private static void help(Options options) {
		HelpFormatter formater = new HelpFormatter();
		formater.printHelp("FastNN", options);
		//System.exit(0);
	}

	public static void main(String[] args) throws IOException, Exception {
		//TODO: Add unit test cases
		//TODO: Perform analysis for publication
		System.err.println("FastNN Version: 0.3.5");
		System.err.println("Java version: " + System.getProperty("java.version"));
		System.err.println("Max Memory: " + Runtime.getRuntime().maxMemory() / (1024*1024) + " mb");
		//System.err.println(Runtime.getRuntime().availableProcessors());
		int nTaxa = 0;
		int nThreads = 1;
		final String fileName;
		final ExecutorService pool;
		final boolean constrained = true;
		boolean original = false;
		final long maxIterations;
		final boolean useMax;
		final double minGrad;
		final int mult;
		final boolean timeMe;
		final double defaultGradMin = -0.0000001;
		final boolean additive;
		//final String[] names;
		//List<NetMaker.AgglomNode> topology;
		//final NetMaker myNM;
		NetMakerOriginal myNMO;
		NetMakerOriginal.NMMode myMode = NetMakerOriginal.NMMode.CANONICAL;
		//NamesAndDistances myND;
		Options options = new Options();
		Option help = new Option("help", "print this message");
		options.addOption(help);
		Option distFile = OptionBuilder.withArgName("file_location")
										.hasArg()
										.withDescription("The distance file in Phyllip format")
										.create("distFile");
		options.addOption(distFile);
		Option threads = OptionBuilder.withArgName("integer")
										 .hasArg()
										 .withDescription("The number of threads to use.  Default: 1")
										 .create("threads");
		options.addOption(threads);
		Option mode = OptionBuilder.withArgName("string")
								   .hasArg()
								   .withDescription("Determines the algorithm mode to run.  The options are: Canonical, Relaxed, Filter, Random_N, Random_NLOGN, Random_LOGN.  Default: Canonical")
								   .create("mode");
		options.addOption(mode);
		Option multiplier = OptionBuilder.withArgName("integer")
				   .hasArg()
				   .withDescription("For the random mode, this gives the constant multiplier that multiplies the search amount.  Default: 5")
				   .create("mult");
		options.addOption(multiplier);
		options.addOption("order", false, "Outputs the circular order only.");
		options.addOption("additive", false, "Performs an additivity check for the relaxed search strategy.");
//		Option maxIter = OptionBuilder.withArgName("integer")
//				   .hasArg()
//				   .withDescription("The maximum number of iterations for split weight estimation.  Setting this to 0 displays the linear ordering only.  Setting it to -1 uses the gradMin only for convergence.  Setting to 1 uses linear least squares.  Default: 10000")
//				   .create("maxIter");
//		options.addOption(maxIter);
//		Option gradMin = OptionBuilder.withArgName("float")
//				   .hasArg()
//				   .withDescription("This determines the convergence criterion.  It is the least value that the gradient can be.  Default: " + defaultGradMin)
//				   .create("gradMin");
//		options.addOption(gradMin);
		//options.addOption("c", false, "Use constrained numerical optimization to calculate the split weights.");
		options.addOption("time", false, "Show timing results.");
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = null;
		try {
			cmd = parser.parse(options, args);
			// Handle options here...
		} catch (ParseException e) {
			System.err.println( "Parsing failed.  Reason: " + e.getMessage() );
			help(options);
			return;
		}
		if (cmd.hasOption("help")) {
			help(options);
			return;
		}
//		if (cmd.hasOption("gradMin")) {
//			String temp = cmd.getOptionValue("gradMin");
//			if (temp == null) {
//				minGrad = defaultGradMin;
//			} else {
//				minGrad = Double.parseDouble(temp);
//				if (minGrad >= 0) {
//					System.err.println("The gradMin needs to be less than 0.");
//					help(options);
//					return;
//				}
//			}
//		} else {
//			minGrad = defaultGradMin;
//		}
//		if (cmd.hasOption("c")) {
//			constrained = true;
//		} else {
//			constrained = false;
//		}
		if (cmd.hasOption("time")) {
			timeMe = true;
		} else {
			timeMe = false;
		}
		if (cmd.hasOption("additive")) {
			additive = true;
		} else {
			additive = false;
		}
		
//		if (cmd.hasOption("maxIter")) {
//			String temp = cmd.getOptionValue("maxIter");
//			if (temp == null) {
//				maxIterations = 10000;
//				useMax = true;
//			} else {
//				maxIterations = Long.parseLong(temp);
//				if (maxIterations <= 0) {
//					useMax = false;
//				} else {
//					useMax = true;
//				}
//			}
//		} else {
//			maxIterations = 10000;
//			useMax = true;
//		}
		if (!cmd.hasOption("distFile")) {
			System.err.println("The program needs a distance file!!");
			help(options);
			return;
		} else {
			fileName = cmd.getOptionValue("distFile");
		}
		if (cmd.hasOption("threads")) {
			String temp = cmd.getOptionValue("threads");
			if (temp == null) {
				nThreads = 1;
			} else {
				nThreads = Integer.parseInt(temp);
			}
		} else {
			nThreads = 1;
		}
		if (cmd.hasOption("mult")) {
			String temp = cmd.getOptionValue("mult");
			if (temp == null) {
				mult = 5;
			} else {
				mult = Integer.parseInt(temp);
			}
		} else {
			mult = 5;
		}
		
		if (cmd.hasOption("mode")) {
			String temp = cmd.getOptionValue("mode");
			myMode = NetMakerOriginal.NMMode.valueOf(temp.toUpperCase(Locale.ENGLISH));
		} else {
			myMode = NetMakerOriginal.NMMode.CANONICAL;
		}
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			try {
				String data = br.readLine();
				data = data.replaceAll("\\s", "");
				nTaxa = Integer.parseInt(data);
			} catch (IOException IOE) {
				throw new RuntimeException(IOE);
			} finally {
				br.close();
			}
		} catch (FileNotFoundException e) {
		    System.err.println("FileNotFound: " + e.getMessage());
		    help(options);
		    //System.err.println(usage);
		    return;
		}
		System.err.println("Calculating a tree for " + nTaxa + " taxa using " + nThreads + " thread(s).");
//		System.err.println("Using maximum iterations: " + maxIterations + " with gradMin: " + minGrad);
		System.err.println("Getting distances from the file: " + fileName);
		//TODO: Reuse distance and name classes instead of doing it all at once.
		double[][] D = null;
		if (nThreads > 1) {
			pool = Executors.newFixedThreadPool(nThreads);
		} else {
			pool = null;
		}
		//TODO: Create 2D array in the first place, but save the 1D array version for the split weights.
		DistancesAndNames danOrg = new DistancesAndNames(fileName, nTaxa);
//		if (myMode.equals(NetMakerOriginal.NMMode.ORIGINAL)) {
//			int total_nodes = 3*nTaxa - 5;
//			D = new double[total_nodes][total_nodes];
//			for (int i = 0; i < nTaxa; i++) {
//				for (int j = 0; j < nTaxa; j++) {
//					D[i+1][j+1] = danOrg.get(i, j);
//				}
//			}
//		} else {
			D = new double[nTaxa][nTaxa];
			for (int i = 0; i < nTaxa; i++) {
				for (int j = 0; j < nTaxa; j++) {
					D[i][j] = danOrg.get(i, j);
				}
			}
//		}
		danOrg = null;
		//myND = null;
		switch (myMode) {
//		case FILTER:
//			//myNM = new NeighborNetFilter(myND.distances, nTaxa, nThreads, pool);
//			myNMO = new NeighborNetFilter(D, nTaxa, nThreads, pool, true);
//			System.err.println("Using the d-filter.");
//			//System.err.println("UNIMPLEMENTED!!!");
//			//return;
//			break;
		case CANONICAL:
			//myNM = new NeighborNetCanonical(myND.distances, nTaxa, nThreads, pool);
			myNMO = new NeighborNetCanonical(D, nTaxa, nThreads, pool);
			System.err.println("Using the canonical implementation.");
			break;
		case RELAXED:
			//myNM = new NeighborNetRelaxed(myND.distances, nTaxa, nThreads, pool);
			myNMO = new NeighborNetLocal(D, nTaxa, nThreads, additive, pool);
			String addString;
			if (additive) {
				addString = " with additivity checking.";
			} else {
				addString = " without additivity checking.";
			}
			System.err.println("Using the relaxed version" + addString);
			break;
//		case ORIGINAL:
//			//TODO: Integrate the original implementation into the NetMaker framework.
//			//NeighborNetOriginal NNO = new NeighborNetOriginal();
//			original = true;
//			myNMO = null;
//			//myNMO = new NeighborNetOriginal(D, nTaxa, nThreads, pool);
//			System.err.println("Using the original implementation.");
//			break;
		case RANDOM_N:
			myNMO = new NeighborNetRandom(D, nTaxa, nThreads, pool, NeighborNetRandom.RandomAmount.N, mult);
			System.err.println("Using the random version with order N search and multiplier " + mult + ".");
			break;
		case RANDOM_NLOGN:
			myNMO = new NeighborNetRandom(D, nTaxa, nThreads, pool, NeighborNetRandom.RandomAmount.NLOGN, mult);
			System.err.println("Using the random version with order NLOGN search and multiplier " + mult + ".");
			break;
		case RANDOM_LOGN:
			myNMO = new NeighborNetRandom(D, nTaxa, nThreads, pool, NeighborNetRandom.RandomAmount.LOGN, mult);
			System.err.println("Using the random version with order LOGN search and multiplier " + mult + ".");
			break;
		default:
			myNMO = new NeighborNetCanonical(D, nTaxa, nThreads, pool);
			System.err.println("Using the canonical implementation.");
//			myNM = new NeighborNetCanonical(myND.distances, nTaxa, nThreads, pool);
//			System.err.println("Using the canonical implementation.");
			break;
		}
		final int[] ordering;
			//topology = null;
		long beginTimer = 0, endTimer, diff;
		if (timeMe) {
			//long beginTimer, endTimer, diff;
//			if (original) {
//				beginTimer = System.nanoTime();
//				ordering = NeighborNetOriginal.runNeighborNet(nTaxa, D);
//				endTimer = System.nanoTime();
//			} else {
				beginTimer = System.nanoTime();
				ordering = myNMO.runNeighborNet();
				endTimer = System.nanoTime();
//			}
			diff = endTimer - beginTimer;
			System.err.println("Got the order in (s): " + diff / 1000000000.0);
			//TODO: Remove the count in final release.
//			if (!original && nThreads == 1) {
//				System.err.println("Total number of cluster pairs examined: " + myNMO.count);
//			}
		} else {
//			if (original) {
//				ordering = NeighborNetOriginal.runNeighborNet(nTaxa, D);
//			} else {
				ordering = myNMO.runNeighborNet();
//			}
		}
		if(cmd.hasOption("order")) {
			System.out.println(Arrays.toString(ordering));
			return;
		}
		if (timeMe) { 
			beginTimer = System.nanoTime();
		}
		int dim = (nTaxa * (nTaxa - 1)) /2;
		//int[] od1 = {0, 1, 2, 3, 4, 5};
		//ordering = od1;
		NonNegativeLeastSquares mySolver = new NonNegativeLeastSquares(dim, dim);
		List<BitSet> splitList = new ArrayList<BitSet>(nTaxa);
		int[] taxa1 = new int[dim];
		int[] taxa2 = new int[dim];
		int count = 0;
		for (int i = 0; i < nTaxa; i++) {
        	BitSet split = new BitSet(nTaxa);
        	for (int j = i + 1; j < nTaxa; j++) {
        		taxa1[count] = i+1;
        		taxa2[count] = j+1;
        		count++;
        		split.set(ordering[j]); //Should this be j+1?
        		splitList.add((BitSet) split.clone());
        	}
        	
        }
//		for (BitSet split : splitList) {
//			System.err.println(split);
//		}
//		System.err.println(Arrays.toString(taxa1));
//		System.err.println(Arrays.toString(taxa2));
		
		
		//double[][] A = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				BitSet mySplit = splitList.get(j);
				if (mySplit.get(taxa1[i]) == mySplit.get(taxa2[i])) {
					mySolver.a[i][j] = 0;
				} else {
					mySolver.a[i][j] = 1;
				}
			}
		}
		final DistancesAndNames dan = new DistancesAndNames(fileName, nTaxa);
		for (int i = 0; i < dan.distances.length; i++) {
			mySolver.b[i] = dan.distances[i];
		}
//		System.err.println("Distances");
//		System.err.println(Arrays.toString(dan.distances));
//		System.err.println("A");
//		for (int i = 0; i < mySolver.M; i++) {
//			System.err.println(Arrays.toString(mySolver.a[i]));
//		}
		//mySolver.a = A;
		//mySolver.b = dan.distances;
		//mySolver.x = new double[dim];
		mySolver.solve();
		//System.err.println(Arrays.toString(mySolver.x));
		double[] x = mySolver.x;
		mySolver = null;
		double optionThreshold = 0.000001;
		List<SplitAndWeight> splits = new ArrayList<SplitAndWeight>(nTaxa);
		int index = 0;
		for (BitSet split : splitList) {
			if (x[index] > optionThreshold) {
				SplitAndWeight saw = new SplitAndWeight();
				saw.split = split;
				saw.weight = x[index];
				splits.add(saw);
			}
			index++;
		}
//        for (int i = 0; i < nTaxa; i++) {
//        	BitSet split = new BitSet(nTaxa);
//        	for (int j = i + 1; j < nTaxa; j++) {
//        		split.set(ordering[j]); //Should this be j+1?
//        		if (x[index] > optionThreshold) {
//        			SplitAndWeight saw = new SplitAndWeight();
//        			saw.split = (BitSet) split.clone();
//        			saw.weight = x[index];
//        			splits.add(saw);
//        		}
//        		index++;
//        	}
//        }
		if (timeMe) {
			endTimer = System.nanoTime();
			diff = endTimer - beginTimer;
			System.err.println("Got the splits and weights in (s): " + diff / 1000000000.0);
			beginTimer = System.nanoTime();
		}
        OutputPrinter.NexusWithSplitsAndDistances(ordering, dan, splits);
		if (timeMe) {
			endTimer = System.nanoTime();
			diff = endTimer - beginTimer;
			System.err.println("Wrote the output in (s): " + diff / 1000000000.0);
		}
		
//		
//		//mySovler.  dan.distances
//			D = null;
//			myNMO = null;
//			//System.out.println(Arrays.toString(ordering));
//			/*Get the splits and weights */
//			System.err.println("Getting the splits and weights...");
//			if (timeMe) { 
//				beginTimer = System.nanoTime();
//			}
//			double optionThreshold = 0.000001;
//			//final double[] weights;
//			//if (nThreads == 1) {
//			if (maxIterations != 0) {
//				final List<SplitAndWeight> splits;
//				final DistancesAndNames dan = new DistancesAndNames(fileName, nTaxa);
//				CircularSplitWeights csw = new CircularSplitWeights();
//				//splits = csw.getWeightedSplits(ordering, dan, "ols", constrained, optionThreshold);
//				splits = csw.getWeightedSplits(ordering, dan, "ols", constrained, optionThreshold, maxIterations, minGrad, useMax);
//				//Splits mySplits = CircularSplitWeights2.getWeightedSplits(ordering, dan, "ols", true, optionThreshold);
//				//} else {
//				//	CircularSplitWeightsMulti csw = new CircularSplitWeightsMulti(nTaxa, nThreads, pool);
//				//	splits = csw.getWeightedSplits(ordering, dan, "ols", constrained, optionThreshold);
//				//}
//				if (timeMe) {
//					endTimer = System.nanoTime();
//					diff = endTimer - beginTimer;
//					System.err.println("Got the splits and weights in (s): " + diff / 1000000000.0);
//					beginTimer = System.nanoTime();
//				}
//				/* Write splits, weights, and data to a nexus file format. */
//				
//				
//				OutputPrinter.NexusWithSplitsAndDistances(ordering, dan, splits);
//				if (timeMe) {
//					endTimer = System.nanoTime();
//					diff = endTimer - beginTimer;
//					System.err.println("Wrote the output in (s): " + diff / 1000000000.0);
//				}
//			} else {
//				System.out.println(Arrays.toString(ordering));
//			}
//			
//		if (pool != null) {
//			pool.shutdown();
//		}
//		System.err.println("Finished!");
		return;
	}
}