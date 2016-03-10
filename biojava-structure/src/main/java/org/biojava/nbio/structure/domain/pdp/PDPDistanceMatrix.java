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

public class PDPDistanceMatrix {
	int[][] dist;
	int nclose;
	int[] iclose ;
	int[] jclose ;

	public PDPDistanceMatrix(){

	}

	public int[][] getDist() {
		return dist;
	}

	public void setDist(int[][] dist) {
		this.dist = dist;
	}

	public int getNclose() {
		return nclose;
	}

	public void setNclose(int nclose) {
		this.nclose = nclose;
	}

	public int[] getIclose() {
		return iclose;
	}

	public void setIclose(int[] iclose) {
		this.iclose = iclose;
	}

	public int[] getJclose() {
		return jclose;
	}

	public void setJclose(int[] jclose) {
		this.jclose = jclose;
	}




}
