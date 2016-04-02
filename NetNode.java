package nnet;

/* A node in the net */

class NetNode implements Comparable<NetNode> {
    int id = 0;  // Mapping to node id in the execution graph
    int distID = 0;  // Mapping to the distance matrix
    int positionID;  // Back mapping to the index in the NetNode array
    NetNode nbr = null; // adjacent node
    NetNode ch1 = null; // first child
    NetNode ch2 = null; // second child
    NetNode next = null; // next in list of active nodes
    NetNode prev = null; // prev in list of active nodes
    //double Rx = 0;  // Row sum
    double Sx = 0;  // Cluster row sum

    public String toString() {
        String str = "[id=" + id;
        str += " distID=" + distID;
        str += " positionID=" + positionID;
        str += " nbr=" + (nbr == null ? "null" : ("" + nbr.id));
        str += " ch1=" + (ch1 == null ? "null" : ("" + ch1.id));
        str += " ch2=" + (ch2 == null ? "null" : ("" + ch2.id));
        str += " prev=" + (prev == null ? "null" : ("" + prev.id));
        str += " next=" + (next == null ? "null" : ("" + next.id));
        //str += " Rx=" + Rx;
        str += " Sx=" + Sx;
        str += "]";
        return str;
    }
    
    public int compareTo(NetNode n) {
    	if (this.Sx < n.Sx) {
    		return -1;
    	} else if (this.Sx > n.Sx) {
    		return 1;
    	} else {
    		return 0;
    	}
    }
    
//    synchronized final void addRx(double adder) {
//    	Rx += adder;
//    }
    
    synchronized final void addSx(double adder) {
    	Sx += adder;
    }
}