package org.biojava.nbio.structure.secstruc;

/** 
 * A Ladder is a set of one or more consecutive bridges of identical type.
 * A Bridge is a Ladder of length one.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public class Ladder {
	
	int from; 	//start of the first strand
	int to;		//end of the first strand
	int lfrom;	//start of the second strand
	int lto;	//end of the second strand
	
	BridgeType btype;
	
	int connectedTo; //another ladder with higher index connected to this
	int connectedFrom; //another ladder with lower index connected to this
	
	public int getFrom() {
		return from;
	}
	public int getTo() {
		return to;
	}
	public int getLfrom() {
		return lfrom;
	}
	public int getLto() {
		return lto;
	}
	public BridgeType getBtype() {
		return btype;
	}
	public int getConnectedTo() {
		return connectedTo;
	}
	public int getConnectedFrom() {
		return connectedFrom;
	}
	
}
