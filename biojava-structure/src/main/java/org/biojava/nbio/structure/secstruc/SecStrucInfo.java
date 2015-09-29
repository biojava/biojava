package org.biojava.nbio.structure.secstruc;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.align.multiple.AbstractScoresCache;

/**
 * Container for the SS information of a single residue.
 * This class is designed to be stored inside an Amino Acid object.
 * It can also contain a back-reference to its parent AA.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucInfo extends AbstractScoresCache {

	/** Secondary strucuture assigned by the PDB author */
	public static final String PDB_AUTHOR_ASSIGNMENT = "PDB_AUTHOR_ASSIGNMENT";
	
	/** Secondary strucuture predicted by the DSSP program */
	public static final String DSSP_ASSIGNMENT = "DSSP_ASSIGNMENT";
	
	private final SecStrucType type;
	private final String assignment;
	private AminoAcid parent;
	
	public SecStrucInfo(AminoAcid parent,String assignment,SecStrucType type){
		super();
		this.type = type;
		this.assignment = assignment;
		this.parent = parent;
	}

	public SecStrucType getType() {
		return type;
	}

	public String getAssignment() {
		return assignment;
	}
	
	public AminoAcid getParent() {
		return parent;
	}
	
}
