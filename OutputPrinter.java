package nnet;

import java.util.List;
import java.util.BitSet;

final public class OutputPrinter {
	
	public static void NexusWithSplitsAndDistances(int[] order, DistancesAndNames dan, List<SplitAndWeight> splits) {
		System.out.println("#nexus");
		System.out.println("");
		PrintTaxa(dan);
		System.out.println("");
		PrintDistances(dan);
		System.out.println("");
		PrintSplits(dan.getNtax(), order, splits);
		System.out.println("");
		PrintAssumptions(dan.getNtax());
		System.out.println("");
	}
	
	private static void PrintTaxa(DistancesAndNames dan) {
		System.out.println("BEGIN Taxa;");
		System.out.println("DIMENSIONS ntax=" + dan.getNtax() + ";");
		System.out.println("TAXLABELS");
		for(int i = 0; i < dan.getNtax(); i++) {
			int j = i + 1;
			String taxaName = dan.names[i];
			System.out.println("[" + j + "]" + " '" + taxaName + "'");
		}
		System.out.println(";");
		System.out.println("END; [Taxa]");
	}
	
	private static void PrintDistances(DistancesAndNames dan) {
		System.out.println("BEGIN Distances;");
		System.out.println("DIMENSIONS ntax=" + dan.getNtax() + ";");
		System.out.println("FORMAT labels=no diagonal triangle=both;");
		System.out.println("MATRIX");
		for(int i = 0; i < dan.getNtax(); i++) {
			for (int j = 0; j < dan.getNtax(); j++) {
				System.out.print(" " + dan.get(i, j));
			}
			System.out.print("\n");
		}
		System.out.println(";");
		System.out.println("END; [Distances]");
	}
	
	private static void PrintSplits(final int ntax, final int[] order, final List<SplitAndWeight> splits) {
		System.out.println("BEGIN Splits;");
		System.out.println("DIMENSIONS ntax=" + ntax + " nsplits=" + splits.size() + ";");
		System.out.println("FORMAT labels=no weights=yes confidences=no intervals=no;");
		System.out.println("PROPERTIES fit=-1.0 cyclic;");
		System.out.print("CYCLE");
		for (int i = 1; i < order.length; i++) {
			int taxa = order[i];
			System.out.print(" " + taxa);
		}
		System.out.print(";\n");
		System.out.println("MATRIX");
		int counter = 1;
		for (SplitAndWeight saw : splits) {
			PrintSplit(counter, ntax, saw);
			counter++;
		}
		System.out.println(";");
		System.out.println("END; [Splits]");
	}
	
	private static void PrintSplit(int i, int ntax, SplitAndWeight saw) {
		BitSet split = saw.split;
		double weight = saw.weight;
		int size = split.cardinality();
		if (ntax - size < size) {
			size = ntax - size;
		}
		System.out.print("[" + i + ", size=" + size + "] \t " + weight + " \t ");
		int taxa = 0;
		for (int j = 0; j < split.cardinality(); j++) {
			taxa = split.nextSetBit(taxa);
			System.out.print(" " + taxa);
			taxa++;
		}
		System.out.print(",\n");
	}
	
	private static void PrintAssumptions(int ntax) {
		System.out.println("BEGIN st_Assumptions;");
		System.out.println("uptodate;");
		System.out.println("disttransform=NeighborNet;");
		System.out.println("splitstransform=EqualAngle;");
		System.out.println("SplitsPostProcess filter=dimension value=" + ntax + ";");
		System.out.println(" exclude  no missing;");
		System.out.println("autolayoutnodelabels;");
		System.out.println("END; [st_Assumptions]");
	}
	
}