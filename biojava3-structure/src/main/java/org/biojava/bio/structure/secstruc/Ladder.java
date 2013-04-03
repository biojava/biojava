package org.biojava.bio.structure.secstruc;

/** Connects two connected strand regions to form larger sheets
 * 
 * @author andreas
 *
 */
public class Ladder {
	int from ;
	int to;
	int lfrom;
	int lto;
	BridgeType btype;
	int connectedTo;
	int connectedFrom;

	public int getFrom() {
		return from;
	}
	public void setFrom(int from) {
		this.from = from;
	}
	public int getTo() {
		return to;
	}
	public void setTo(int to) {
		this.to = to;
	}
	public int getLfrom() {
		return lfrom;
	}
	public void setLfrom(int lfrom) {
		this.lfrom = lfrom;
	}
	public int getLto() {
		return lto;
	}
	public void setLto(int lto) {
		this.lto = lto;
	}
	public BridgeType getBtype() {
		return btype;
	}
	public void setBtype(BridgeType btype) {
		this.btype = btype;
	}
	public int getConnectedTo() {
		return connectedTo;
	}
	public void setConnectedTo(int connectedTo) {
		this.connectedTo = connectedTo;
	}
	public int getConnectedFrom() {
		return connectedFrom;
	}
	public void setConnectedFrom(int connectedFrom) {
		this.connectedFrom = connectedFrom;
	}
	@Override
	public String toString() {
		return "Ladder [from=" + from + ", to=" + to + ", lfrom=" + lfrom
				+ ", lto=" + lto + ", btype=" + btype + ", connectedTo="
				+ connectedTo + ", connectedFrom=" + connectedFrom + "]";
	}
	

}
