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
 * Created on Jun 17, 2010
 * Author: ap3 
 *
 */

package org.biojava.bio.structure;

import java.io.Serializable;
import java.io.StringWriter;

/** Everything that is needed to uniquely describe a residue position
 * 
 * @author Andreas Prlic
 *
 */
public class ResidueNumber implements Serializable, Comparable<ResidueNumber>
{

	private static final long serialVersionUID = 1773011704758536083L;
	private String chainId;
	private Character insCode;
	private Integer seqNum;

	public ResidueNumber() {
	}

	public ResidueNumber(String chainId, Integer residueNumber, Character insCode) {
		this.chainId = chainId;
		this.seqNum = residueNumber;
		this.insCode = insCode;
	}

	public String getChainId()
	{
		return chainId;
	}
	public void setChainId(String chainId)
	{
		this.chainId = chainId;
	}
	public Character getInsCode()
	{
		return insCode;
	}
	public void setInsCode(Character insCode)
	{
		this.insCode = insCode;
	}
	public Integer getSeqNum()
	{
		return seqNum;
	}
	public void setSeqNum(Integer seqNum)
	{
		this.seqNum = seqNum;
	}

	/**
	 * Require chainId, insCode, and seqNum to be equal
	 * @param obj
	 * @return
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ResidueNumber))
			return false;

		if (obj==this)
			return true;

		ResidueNumber anNumber = (ResidueNumber) obj;



		if (insCode!=null) {
			if (insCode != anNumber.getInsCode()) 
				return false;
		} else {
			if (anNumber.getInsCode()!=null)
				return false;
		}
		if (seqNum!=null) {
			if (!seqNum.equals(anNumber.getSeqNum()))
				return false;
		} else {
			if (anNumber.getSeqNum()!=null)
				return false;
		}

		if ( chainId != null) {
			if ( ! chainId.equals(anNumber.getChainId()))
				return false;
		} else {
			if ( anNumber.getChainId() != null)
				return false;
		}

		return true;
	}

	@Override
	public int hashCode() {
		int result = 17;
		result = 31 * result + (chainId != null ? chainId.hashCode() : 0);
		result = 31 * result + (seqNum != null ? seqNum.hashCode() : 0);
		result = 31 * result + (insCode==null ? 0 : insCode.hashCode());
		return result;
	}

	/**
	 * @return The residue number and insertion code as a string, eg "74A"
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {

		StringWriter writer = new StringWriter();
		//	   if ( chainId != null){
			//		   writer.append(chainId);
			//		   writer.append(":");
			//	   }
		writer.append(seqNum+"");
		if (  insCode != null && 
				( insCode != ' '))
			writer.append(insCode);

		return writer.toString();
	}

	/**
	 * @return The chain, number, and insertion code as a string, eg "B  74A" or "A    1 "
	 */
	public String toPDB() {
		String insCodeS ;
		if ( insCode != null)
			insCodeS = insCode+"";
		else insCodeS = " ";
		return String.format("%s%4d%-2s", chainId, seqNum, insCodeS);
	}

	
	/** Convert a string representation of a residue number to a residue number object.
	 * The string representation can be a integer followed by a character.
	 * 
	 * @param pdb_code
	 * @return a ResidueNumber object
	 */
	public static ResidueNumber fromString(String pdb_code) {
		ResidueNumber residueNumber = new ResidueNumber();
		Integer resNum = null;
		String icode = null;		
		
		try {
			resNum = Integer.parseInt(pdb_code);
		} catch ( NumberFormatException e){
			// there is an insertion code..
		
			// Split at any position that's either:
			// preceded by a digit and followed by a non-digit, or
			// preceded by a non-digit and followed by a digit.
			String[] spl = pdb_code.split("(?<=\\d)(?=\\D)|(?<=\\D)(?=\\d)");
			if ( spl.length == 2){
				resNum = Integer.parseInt(spl[0]);
				icode = spl[1];
			}
			
		}
				
		residueNumber.setSeqNum(resNum);
		if ( icode == null)
			residueNumber.setInsCode(null);
		else if ( icode.length() > 0)
			residueNumber.setInsCode(icode.charAt(0));
		return residueNumber;
	}

	@Override
	public int compareTo(ResidueNumber other) {
	
		// compare chain ID
		if ( this.getChainId() != null && other.getChainId() != null ) {
			int chainComp = this.getChainId().compareTo(other.getChainId()); 
			if (  chainComp != 0 )
				return chainComp;
			
		}
			
		// if equal sequence position, compare insCode
		if ( this.getSeqNum() == other.getSeqNum()) {
			if (this.getInsCode() != null && other.getInsCode() != null) {
				return this.getInsCode().compareTo(other.getInsCode());
			}
			if ( this.getInsCode() == null && other.getInsCode() != null)
				return -1;
			else return 1;
		}
		return getSeqNum().compareTo(other.getSeqNum());
		
	}
}
