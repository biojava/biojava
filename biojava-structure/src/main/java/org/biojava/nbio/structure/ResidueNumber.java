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

package org.biojava.nbio.structure;

import java.io.Serializable;
import java.io.StringWriter;

/** 
 * Everything that is needed to uniquely describe a residue position
 *
 * @author Andreas Prlic
 *
 */
public class ResidueNumber implements Serializable, Comparable<ResidueNumber>
{

	private static final long serialVersionUID = 1773011704758536083L;
	private String chainName;
	private Character insCode;
	private Integer seqNum;

	public ResidueNumber() {
	}

	public ResidueNumber(ResidueNumber o) {
		this.chainName = o.chainName;
		this.insCode = o.insCode;
		this.seqNum = o.seqNum;
	}

	public ResidueNumber(String chainName, Integer residueNumber, Character insCode) {
		this.chainName = chainName;
		this.seqNum = residueNumber;
		this.insCode = insCode;
	}

	public String getChainName()
	{
		return chainName;
	}
	public void setChainName(String chainName)
	{
		this.chainName = chainName;
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
		if (chainName == null) {
			if (other.chainName != null)
				return false;
		} else if (!chainName.equals(other.chainName))
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
	
	/**
	 * Check if the seqNum and insertion code are equivalent,
	 * ignoring the chain
	 * @param obj
	 * @return
	 */
	public boolean equalsPositional(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ResidueNumber other = (ResidueNumber) obj;
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
		result = prime * result + ((chainName == null) ? 0 : chainName.hashCode());
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
		//	   if ( chainName != null){
		//		   writer.append(chainName);
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
		return String.format("%s%4d%-2s", chainName, seqNum, insCodeS);
	}


	/** Convert a string representation of a residue number to a residue number object.
	 * The string representation can be a integer followed by a character.
	 *
	 * @param pdb_code
	 * @return a ResidueNumber object, or null if the input was null
	 */
	public static ResidueNumber fromString(String pdb_code) {
		if(pdb_code == null)
			return null;

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


	/**
	 * Compare residue numbers by chain, sequence number, and insertion code
	 */
	@Override
	public int compareTo(ResidueNumber other) {

		// chain id
		if (chainName != null && other.chainName != null) {
			if (!chainName.equals(other.chainName)) return chainName.compareTo(other.chainName);
		}
		if (chainName != null && other.chainName == null) {
			return 1;
		} else if (chainName == null && other.chainName != null) {
			return -1;
		}

		return compareToPositional(other);
	}

	/**
	 * Compare residue numbers by sequence number and insertion code,
	 * ignoring the chain
	 * @param other
	 * @return
	 */
	public int compareToPositional(ResidueNumber other) {
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
		final String chain = chainName==null? "" : chainName;
		return chain + "_" + toString();
	}

}
