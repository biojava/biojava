/**
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
 * Created on 2013-06-13
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.nbio.structure.rcsb;

import java.util.ArrayList;
import java.util.List;

/**
 * Corresponds to the wrapper element "ligandInfo" in an RCSB {@code ligandInfo} XML file.
 *
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 *
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBLigands {

	private String pdbId;

	private List<RCSBLigand> ligands;

	public RCSBLigands() {
		ligands = new ArrayList<RCSBLigand>();
	}

	public void addLigand(RCSBLigand ligand) {
		ligands.add(ligand);
	}

	public String getPdbId() {
		return pdbId;
	}

	public List<RCSBLigand> getLigands() {
		return ligands;
	}

	void setPdbId(String pdbId) {
		this.pdbId = pdbId;
	}

}
