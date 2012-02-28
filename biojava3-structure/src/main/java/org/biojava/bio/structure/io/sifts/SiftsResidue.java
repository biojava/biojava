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
package org.biojava.bio.structure.io.sifts;

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
				+ pdbResName + ", chainId=" + chainId + ", uniProtResName="
				+ uniProtResName + ", uniProtPos=" + uniProtPos
				+ ", naturalPos=" + naturalPos + ", seqResName=" + seqResName
				+ ", pdbId=" + pdbId + ", uniProtAccessionId="
				+ uniProtAccessionId + ", notObserved=" + notObserved + "]";
	}

	

	
}
