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
package org.biojava.nbio.structure.domain.pdp;



public class CutValues {
	public double s_min;
	public int site2;
	public boolean first_cut;

	public double AD;

	public CutValues(){
		s_min = 100;
		site2 = 0;
		first_cut = true;
	}

	@Override
	public String toString() {
		return "CutValues [s_min=" + s_min + ", site2=" + site2 +
		", AD=" + AD
		+ "]";
	}



}
