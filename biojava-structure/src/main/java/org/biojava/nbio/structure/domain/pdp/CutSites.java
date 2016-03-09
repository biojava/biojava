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

import java.util.Arrays;


public class CutSites {
	int ncuts;
	int[] cut_sites;

	public CutSites(){
		ncuts = 0;

		cut_sites = new int[PDPParameters.MAX_CUTS];
	}

	@Override
	public String toString() {
		return "CutSites [ncuts=" + ncuts + ", cut_sites="
				+ Arrays.toString(cut_sites) + "]";
	}
	public int getNcuts() {
		return ncuts;
	}
	public void setNcuts(int ncuts) {
		this.ncuts = ncuts;
	}
	public int[] getCut_sites() {
		return cut_sites;
	}
	public void setCut_sites(int[] cut_sites) {
		this.cut_sites = cut_sites;
	}



}
