/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.secstruc;

/** 
 * Connects two connected strand regions to form larger sheets
 * 
 * @author Andreas Prlic
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
