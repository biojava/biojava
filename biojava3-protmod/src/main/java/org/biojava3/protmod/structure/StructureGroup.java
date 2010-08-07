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
 * Created on Aug 2, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.structure;

import java.io.Serializable;

import org.biojava.bio.structure.PDBResidueNumber;

import org.biojava3.protmod.ComponentType;

/**
 * Information of a group (residue or ligand) involved in a modification. 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class StructureGroup
implements Serializable {
	private static final long serialVersionUID = -5648208521422258019L;
	
	private final PDBResidueNumber resNum;
	private final String pdbCode;
	private final ComponentType type;
	
	public StructureGroup(final PDBResidueNumber resNum,
			final String pdbCode, final ComponentType type) {
		this.resNum = resNum;
		this.pdbCode = pdbCode;
		this.type = type;
	}

	public PDBResidueNumber getPDBResidueNumber() {
		return resNum;
	}
	
	public String getChainId() {
		return resNum.getChainId();
	}
	
	public int getResidueNumber() {
		return resNum.getResidueNumber();
	}
	
	public String getInsCode() {
		return resNum.getInsCode();
	}

	public String getPDBCode() {
		return pdbCode;
	}

	public ComponentType getType() {
		return type;
	}
	
	public boolean isAminoAcid() {
		return type == ComponentType.AMINOACID;
	}
	
	public boolean equals(Object obj) {
		if (obj == this)
			return false;
		
		if (!(obj instanceof StructureGroup))
			return false;
		
		StructureGroup aGroup = (StructureGroup) obj;
		if (!resNum.equals(aGroup.resNum))
			return false;
		
		if (type != aGroup.type)
			return false;
		
		return true;
	}
	
	public int hashCode() {
		int result = 17;
		result = result * 31 + resNum.hashCode();
		result = result * 31 + type.hashCode();
		return result;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(pdbCode);
		sb.append('\t');
		sb.append(resNum.getChainId());
		sb.append('\t');
		sb.append(resNum.getResidueNumber());
		if (resNum.getInsCode() != null)
			sb.append(resNum.getInsCode());
		sb.append('\t');
		return sb.toString();
	}
}
