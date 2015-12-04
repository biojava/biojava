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

import org.biojava.nbio.structure.Group;

/**
 * Container for the secondary structure information of a single residue. This
 * class is designed to be stored inside an Amino Acid object. It can also
 * contain a back-reference to its parent AA.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucInfo {

	/** Secondary strucuture assigned by the PDB author */
	public static final String PDB_AUTHOR_ASSIGNMENT = "PDB_AUTHOR_ASSIGNMENT";

	/** Secondary strucuture parsed from a DSSP output file */
	public static final String DSSP_ASSIGNMENT = "DSSP_ASSIGNMENT";

	/** Secondary strucuture calculated and assigned by DSSP of BioJava */
	public static final String BIOJAVA_ASSIGNMENT = "BIOJAVA_ASSIGNMENT";

	protected SecStrucType type;
	protected String assignment;
	protected Group parent;

	public SecStrucInfo(Group g, String ass, SecStrucType t) {
		type = t;
		assignment = ass;
		parent = g;
	}

	public SecStrucType getType() {
		return type;
	}

	public void setType(SecStrucType t) {
		type = t;
	}

	public String getAssignment() {
		return assignment;
	}

	public Group getGroup() {
		return parent;
	}

	@Override
	public String toString() {
		return assignment + ": " + type;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SecStrucInfo))
			return false;
		else {
			SecStrucInfo ss = (SecStrucInfo) o;
			if (type == ss.type)
				return true;
			else
				return false;
		}
	}

}
