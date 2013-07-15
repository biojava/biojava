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
 * Created on 2012-11-20
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.rcsb;

import java.util.ArrayList;
import java.util.List;

/**
 * Corresponds to the wrapper element in an RCSB {@code describeMol} XML file.
 * 
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 * 
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBDescription {

	private String pdbId;

	private List<RCSBPolymer> polymers;

	public RCSBDescription() {
		polymers = new ArrayList<RCSBPolymer>();
	}

	public void addPolymer(RCSBPolymer polymer) {
		polymers.add(polymer);
	}

	public String getPdbId() {
		return pdbId;
	}

	public List<RCSBPolymer> getPolymers() {
		return polymers;
	}

	void setPdbId(String pdbId) {
		this.pdbId = pdbId;
	}

}
