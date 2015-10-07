package org.biojava.nbio.structure.secstruc;

import org.biojava.nbio.structure.Group;

/**
 * Container for the secondary structure information of a single residue.
 * This class is designed to be stored inside an Amino Acid object.
 * It can also contain a back-reference to its parent AA.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucInfo {

	/** Secondary strucuture assigned by the PDB author */
	public static final String PDB_AUTHOR_ASSIGNMENT = "PDB_AUTHOR_ASSIGNMENT";
	
	/** Secondary strucuture parsed from a DSSP output file */
	public static final String DSSP_FILE_ASSIGNMENT = "DSSP_ASSIGNMENT";
	
	/** Secondary strucuture calculated and assigned by DSSP of BioJava */
	public static final String BIOJAVA_ASSIGNMENT = "BIOJAVA_ASSIGNMENT";
	
	protected SecStrucType type;
	protected String assignment;
	protected Group parent;
	
	public SecStrucInfo(Group g, String ass, SecStrucType t){
		type = t;
		assignment = ass;
		parent = g;
	}

	public SecStrucType getType() {
		return type;
	}
	
	public void setType(SecStrucType t){
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
	public boolean equals(Object o){
		if (!(o instanceof SecStrucInfo)) return false;
		else {
			SecStrucInfo ss = (SecStrucInfo) o;
			if (type == ss.type) return true;
			else return false;
		}
	}
	
}
