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
 * Corresponds to a macromolecule in an RCSB {@code describeMol} XML file.
 * 
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 * 
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBMacromolecule {

	private List<String> accessions;
	private String name;

	public RCSBMacromolecule() {
		accessions = new ArrayList<String>();
	}

	public List<String> getAccessions() {
		return accessions;
	}

	public String getName() {
		return name;
	}

	void addAccession(String e) {
		accessions.add(e);
	}

	void setName(String name) {
		this.name = name;
	}

}
