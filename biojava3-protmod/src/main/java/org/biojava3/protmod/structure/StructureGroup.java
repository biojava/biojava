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

import org.biojava.bio.structure.ResidueNumber;

/**
 * Information of a group (residue or ligand) involved in a modification. 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class StructureGroup
implements Comparable<StructureGroup> {
	
	private  ResidueNumber resNum;
	private  String pdbName;
        private Boolean isAminoAcid;
	
	public StructureGroup(){
		resNum = new ResidueNumber();
	}

	public StructureGroup( ResidueNumber resNum,
			 String pdbName, boolean isAminoAcid) {
		this.resNum = resNum;
		this.pdbName = pdbName;
        this.isAminoAcid = isAminoAcid;
	}

//	public StructureGroup( ResidueNumber resNum,
//			 String pdbName) {
//		this(resNum, pdbName, true);
//	}

	public ResidueNumber getPDBResidueNumber() {
		return resNum;
	}
	
	public void setPDBResidueNumber(ResidueNumber resNum) {
		this.resNum = resNum;
	}
	public String getChainId() {
		return resNum.getChainId();
	}
	
	public void setChainId(String chainId){
		if ( resNum == null)
			resNum = new ResidueNumber();
		resNum.setChainId(chainId);
	}
	
	public int getResidueNumber() {
		return resNum.getSeqNum();
	}
	
	public void setResidueNumber(int seqNr){
		if ( resNum == null)
			resNum = new ResidueNumber();
		resNum.setSeqNum(seqNr);
	}
	
	public Character getInsCode() {
		return resNum.getInsCode();
	}
	
	public void setInsCode(Character c){
		if ( resNum == null)
			resNum = new ResidueNumber();
		resNum.setInsCode(c);
	}

	public String getPDBName() {
		return pdbName;
	}
	
	public void setPDBName(String pdbName){
		this.pdbName = pdbName;
		
	}
	
	public void setIsAminoAcid(boolean isAminoAcid) {
		this.isAminoAcid = isAminoAcid;
	}

	public boolean isAminoAcid() {
		return isAminoAcid;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == this)
			return true;
		
		if (!(obj instanceof StructureGroup))
			return false;
		
		StructureGroup aGroup = (StructureGroup) obj;
		if (!resNum.equals(aGroup.resNum))
			return false;
		
		return true;
	}

	@Override
	public int hashCode() {
		int result = 17;
		if ( resNum != null)
			result = result * 31 + resNum.hashCode();
		return result;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(pdbName);
		sb.append('\t');
		sb.append(resNum.getChainId());
		sb.append('\t');
		sb.append(resNum.getSeqNum());
		if (resNum.getInsCode() != null)
			sb.append(resNum.getInsCode());
		sb.append('\t');
		return sb.toString();
	}

	@Override
	public int compareTo(StructureGroup aGroup) {
		int result = getChainId().compareTo(aGroup.getChainId());
		if (result != 0)
			return result;
		result = getResidueNumber()-aGroup.getResidueNumber();
		if (result != 0)
			return result;
		if (getInsCode()==null) {
			if (aGroup.getInsCode()!=null)
				return -1;
		} else {
			if (aGroup.getInsCode()==null)
				return 1;
			else {
				result = getInsCode().compareTo(aGroup.getInsCode());
				if (result != 0)
					return result;
			}
		}
		return 0;
	}
}
