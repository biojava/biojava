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

	

	

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ResidueNumber other = (ResidueNumber) obj;
		if (chainId == null) {
			if (other.chainId != null)
				return false;
		} else if (!chainId.equals(other.chainId))
			return false;
		if (insCode == null) {
			if (other.insCode != null)
				return false;
		} else if (!insCode.equals(other.insCode))
			return false;
		if (seqNum == null) {
			if (other.seqNum != null)
				return false;
		} else if (!seqNum.equals(other.seqNum))
			return false;
		
		return true;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chainId == null) ? 0 : chainId.hashCode());
		result = prime * result + ((insCode == null) ? 0 : insCode.hashCode());
		result = prime * result + ((seqNum == null) ? 0 : seqNum.hashCode());
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
        writer.append(String.valueOf(seqNum));
		if (  insCode != null && ( insCode != ' '))
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


	public int compareTo(ResidueNumber other) {
		
		// chain id
		if (chainId != null && other.chainId != null) {
			if (!chainId.equals(other.chainId)) return chainId.compareTo(other.chainId);
		}
		if (chainId != null && other.chainId == null) {
			return 1;
		} else if (chainId == null && other.chainId != null) {
			return -1;
		}
		
		// sequence number
		if (seqNum != null && other.seqNum != null) {
			if (!seqNum.equals(other.seqNum)) return seqNum.compareTo(other.seqNum);
		}
		if (seqNum != null && other.seqNum == null) {
			return 1;
		} else if (seqNum == null && other.seqNum != null) {
			return -1;
		}
		
		// insertion code
		if (insCode != null && other.insCode != null) {
			if (!insCode.equals(other.insCode)) return insCode.compareTo(other.insCode);
		}
		if (insCode != null && other.insCode == null) {
			return 1;
		} else if (insCode == null && other.insCode != null) {
			return -1;
		}
		
		return 0;
	}

	public String printFull() {
		final String chain = chainId==null? "" : chainId;
		return chain + toString();
	}
	
}
