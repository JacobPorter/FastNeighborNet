Fast Neighbor-Net is a Java 7 program for computing a phylogenetic split network with the Neighbor-Net approach.

It takes as input a Phylip distance file of taxa.

It produces as output a Nexus format file that can be rendered in SplitsTree4.

Error messages, etc. are printed to stderr and the Nexus file is output to stdout.

To visualize the network, open the file in SplitsTree4 and then choose the menu option Draw -> EqualAngle -> Apply.

The following is an example usage.  It gives the heap 5g, runs the relaxed search strategy with 2 threads and saves the output.
java -server -Xmx5g -jar FastNN.v0.3.5.jar -distFile <path_to_file> -mode relaxed -threads 2 1> <path_to_output>

For large inputs, the Java heap size may need to be increased with the -Xmx option, and it is recommended that java be run in server mode.

The following is the help file for the software.

usage: FastNN
 -additive                   Performs an additivity check for the relaxed
                             search strategy.
 -distFile <file_location>   The distance file in Phyllip format
 -help                       print this message
 -mode <string>              Determines the algorithm mode to run.  The
                             options are: Canonical, Relaxed, 
                             Random_N, Random_NLOGN, Random_LOGN.
                             Default: Canonical
 -mult <integer>             For the random mode, this gives the constant
                             multiplier that multiplies the search amount.
                             Default: 5
 -order                      Outputs the circular order only.
 -threads <integer>          The number of threads to use.  Default: 1
 -time                       Show timing results.

The GNU license applies.