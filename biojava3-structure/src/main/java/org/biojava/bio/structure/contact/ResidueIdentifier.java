package org.biojava.bio.structure.contact;

/**
 * A class used only within the contact package to be able to compare 
 * contacts based on residue numbers and independently from chain identifiers
 * 
 * @author duarte_j
 */
class ResidueIdentifier {

	private int seqNum;
	private Character insCode;
	
	public ResidueIdentifier(int seqNum, Character insCode) {
		this.seqNum = seqNum;
		this.insCode = insCode;
	}
	
	public int getSeqNum() {
		return seqNum;
	}

	public void setSeqNum(int seqNum) {
		this.seqNum = seqNum;
	}

	public Character getInsCode() {
		return insCode;
	}

	public void setInsCode(char insCode) {
		this.insCode = insCode;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + seqNum;
		result = prime * result + (insCode==null?0:insCode);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ResidueIdentifier other = (ResidueIdentifier) obj;
		if (this.seqNum!=other.seqNum) 
			return false;
		if (this.insCode==null && other.insCode!=null)
			return false;
		if (this.insCode!=null && other.insCode==null)
			return false;
		if (this.insCode==null && other.insCode==null)
			return true;
		if (this.insCode!=other.insCode)
			return false;

		return true;
	}

	@Override
	public String toString() {
		return ""+ seqNum + (insCode==null?"":insCode);
	}

	
}
