/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

import org.biojava.bio.structure.ResidueNumber;

/**
 * Small bean to hold information about a single residue in the PDB
 * if we don't want to create a full {@link org.biojava.bio.structure.Group Group} object.
 * @author Spencer Bliven
 * @deprecated replaced by {@link ResidueNumber}
 */
public class PDBResidue {
	private String residueCode;
	private String chain;
	private String aaName; //3-letter code
	
	/**
	 * @param residueCode The {@link org.biojava.bio.structure.Group#getPDBCode() residue code}
	 *  for this residue (residue number + insertion code)
	 * @param chain 1-letter String giving the chain
	 * @param aaName 3-letter String giving the amino acid
	 */
	public PDBResidue(String residueCode, String chain, String aaName) {
		this.residueCode = residueCode;
		this.chain = chain;
		this.aaName = aaName;
	}
	
	/**
	 * @param residueCode The {@link org.biojava.bio.structure.Group#getPDBCode() residue code}
	 *  for this residue (residue number + insertion code)
	 * @param chain 1-letter String giving the chain
	 */
	public PDBResidue(String residueCode, String chain) {
		this(residueCode,chain,null);
	}

	/**
	 * @return The {@link org.biojava.bio.structure.Group#getPDBCode() residue code}
	 *  for this residue (residue number + insertion code)
	 */
	public String getResidueCode() {
		return residueCode;
	}

	/**
	 * @param residueCode the residueCode to set
	 */
	public void setResidueCode(String residueCode) {
		this.residueCode = residueCode;
	}

	/**
	 * @return the chain
	 */
	public String getChain() {
		return chain;
	}

	/**
	 * @param chain the chain to set
	 */
	public void setChain(String chain) {
		this.chain = chain;
	}

	/**
	 * @return the 3-letter amino acid code, or null if none is set
	 */
	public String getAaName() {
		return aaName;
	}

	/**
	 * @param aaName the 3-letter amino acid code
	 */
	public void setAaName(String aaName) {
		this.aaName = aaName;
	}

	/**
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String str = "";
		if(aaName != null) {
			str+=aaName+".";
		}
		return String.format("%s%s.%s", str,residueCode,chain );
	}
	
	
	
	
	

}
