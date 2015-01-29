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

/**
 * Corresponds to a taxonomy in a {@code describeMol} XML file.
 * 
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 * 
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBTaxonomy {

	private final int id;
	private final String name;

	public RCSBTaxonomy(String name, int id) {
		this.name = name;
		this.id = id;
	}

	public int getId() {
		return id;
	}

	public String getName() {
		return name;
	}

}
