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
 * A Ladder is a set of one or more consecutive bridges of identical type. A
 * Bridge is a Ladder of length one.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public class Ladder {

	int from; // start of the first strand
	int to; // end of the first strand
	int lfrom; // start of the second strand
	int lto; // end of the second strand

	BridgeType btype;

	int connectedTo; // another ladder with higher index connected to this
	int connectedFrom; // another ladder with lower index connected to this

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
