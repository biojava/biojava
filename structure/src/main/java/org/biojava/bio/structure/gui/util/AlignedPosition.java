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
 * created at May 28, 2008
 */
package org.biojava.bio.structure.gui.util;

public class AlignedPosition {
	int pos1;
	int pos2;
	int equivalent;

	/** flag if this position is equivalent
	 * 
	 */
	public static final int EQUIVALENT  = 1;
	
	/** they can be shown in the same column (for a compact display)
	 * , but they are not structurally equivalent
	 * 
	 */
	public static final int NOT_ALIGNED = 0;
	
	public AlignedPosition(){
		pos1 = -1;
		pos2 = -1;
		equivalent = NOT_ALIGNED;
	}
	
	public int getPos(int position){
		if (position == 1)
			return pos1;
		else 
			return pos2;
	}
	public void setPos(int position, int pos){
		if (position == 1)
			pos1 = pos;
		else 
			pos2 = pos;
	}
	
	

	
	public  String toString(){
		String t = " AlignedPosition pos1: " + pos1 + " pos2: "+ pos2 ;
		if ( equivalent == EQUIVALENT)
			t+= " EQR";
		
		return t;
	}
	
	public int getPos1() {
		return pos1;
	}

	public void setPos1(int pos1) {
		this.pos1 = pos1;
	}

	public int getPos2() {
		return pos2;
	}

	public void setPos2(int pos2) {
		this.pos2 = pos2;
	}

	public int getEquivalent() {
		return equivalent;
	}

	public void setEquivalent(int equivalent) {
		this.equivalent = equivalent;
	}
	
	
}
