/**
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
 * Created on Feb 22, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.io.sifts;

import java.io.Serializable;

public class SiftsResidue implements Serializable{

	/**
	 *
	 */
	private static final long serialVersionUID = 3425769737629800828L;
	String pdbResNum;
	String pdbResName;
	String chainId;
	String uniProtResName;
	Integer uniProtPos;
	Integer naturalPos;
	String seqResName;
	String pdbId;
	String uniProtAccessionId;
	Boolean notObserved;

	public String getPdbResNum() {
		return pdbResNum;
	}
	public void setPdbResNum(String pdbResNum) {
		this.pdbResNum = pdbResNum;
	}
	public String getPdbResName() {
		return pdbResName;
	}
	public void setPdbResName(String pdbResName) {
		this.pdbResName = pdbResName;
	}
	public String getChainId() {
		return chainId;
	}
	public void setChainId(String chainId) {
		this.chainId = chainId;
	}
	public String getUniProtResName() {
		return uniProtResName;
	}
	public void setUniProtResName(String uniProtResName) {
		this.uniProtResName = uniProtResName;
	}
	public Integer getUniProtPos() {
		return uniProtPos;
	}
	public void setUniProtPos(Integer uniProtPos) {
		this.uniProtPos = uniProtPos;
	}

	public void setPdbId(String dbAccessionId) {
		this.pdbId = dbAccessionId;
	}
	public String getPdbId(){
		return pdbId;
	}
	public void setUniProtAccessionId(String dbAccessionId) {
		this.uniProtAccessionId = dbAccessionId;
	}
	public String getUniProtAccessionId(){
		return uniProtAccessionId;
	}
	public Integer getNaturalPos() {
		return naturalPos;
	}
	public void setNaturalPos(Integer naturalPos) {
		this.naturalPos = naturalPos;
	}
	public Boolean getNotObserved() {
		return notObserved;
	}
	public void setNotObserved(Boolean notObserved) {
		this.notObserved = notObserved;
	}
	public String getSeqResName() {
		return seqResName;
	}
	public void setSeqResName(String seqResName) {
		this.seqResName = seqResName;
	}
	@Override
	public String toString() {
		return "SiftsResidue [pdbResNum=" + pdbResNum + ", pdbResName="
				+ pdbResName + ", chainName=" + chainId + ", uniProtResName="
				+ uniProtResName + ", uniProtPos=" + uniProtPos
				+ ", naturalPos=" + naturalPos + ", seqResName=" + seqResName
				+ ", pdbId=" + pdbId + ", uniProtAccessionId="
				+ uniProtAccessionId + ", notObserved=" + notObserved + "]";
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chainId == null) ? 0 : chainId.hashCode());
		result = prime * result + ((naturalPos == null) ? 0 : naturalPos.hashCode());
		result = prime * result + ((notObserved == null) ? 0 : notObserved.hashCode());
		result = prime * result + ((pdbId == null) ? 0 : pdbId.hashCode());
		result = prime * result + ((pdbResName == null) ? 0 : pdbResName.hashCode());
		result = prime * result + ((pdbResNum == null) ? 0 : pdbResNum.hashCode());
		result = prime * result + ((seqResName == null) ? 0 : seqResName.hashCode());
		result = prime * result + ((uniProtAccessionId == null) ? 0 : uniProtAccessionId.hashCode());
		result = prime * result + ((uniProtPos == null) ? 0 : uniProtPos.hashCode());
		result = prime * result + ((uniProtResName == null) ? 0 : uniProtResName.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		SiftsResidue other = (SiftsResidue) obj;
		if (chainId == null) {
			if (other.chainId != null) return false;
		} else if (!chainId.equals(other.chainId)) return false;
		if (naturalPos == null) {
			if (other.naturalPos != null) return false;
		} else if (!naturalPos.equals(other.naturalPos)) return false;
		if (notObserved == null) {
			if (other.notObserved != null) return false;
		} else if (!notObserved.equals(other.notObserved)) return false;
		if (pdbId == null) {
			if (other.pdbId != null) return false;
		} else if (!pdbId.equals(other.pdbId)) return false;
		if (pdbResName == null) {
			if (other.pdbResName != null) return false;
		} else if (!pdbResName.equals(other.pdbResName)) return false;
		if (pdbResNum == null) {
			if (other.pdbResNum != null) return false;
		} else if (!pdbResNum.equals(other.pdbResNum)) return false;
		if (seqResName == null) {
			if (other.seqResName != null) return false;
		} else if (!seqResName.equals(other.seqResName)) return false;
		if (uniProtAccessionId == null) {
			if (other.uniProtAccessionId != null) return false;
		} else if (!uniProtAccessionId.equals(other.uniProtAccessionId)) return false;
		if (uniProtPos == null) {
			if (other.uniProtPos != null) return false;
		} else if (!uniProtPos.equals(other.uniProtPos)) return false;
		if (uniProtResName == null) {
			if (other.uniProtResName != null) return false;
		} else if (!uniProtResName.equals(other.uniProtResName)) return false;
		return true;
	}

}
