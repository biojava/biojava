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
 */
package org.biojava.nbio.structure;

import java.io.Serializable;

/**
 * Created by andreas on 9/11/15.
 */
public class SeqMisMatchImpl implements SeqMisMatch, Serializable{

	private static final long serialVersionUID = -3699285122925652562L;

	Integer seqNum;
	String origGroup;
	String pdbGroup;
	String details;
	String uniProtId;
	String insCode;
	String pdbResNum;

	@Override
	public Integer getSeqNum() {
		return seqNum;
	}

	@Override
	public void setSeqNum(Integer seqNum) {
		this.seqNum = seqNum;
	}

	@Override
	public String getOrigGroup() {
		return origGroup;
	}

	@Override
	public void setOrigGroup(String origGroup) {
		this.origGroup = origGroup;
	}

	@Override
	public String getPdbGroup() {
		return pdbGroup;
	}

	@Override
	public void setPdbGroup(String pdbGroup) {
		this.pdbGroup = pdbGroup;
	}

	@Override
	public String getDetails() {
		return details;
	}

	@Override
	public void setDetails(String details) {
		this.details = details;
	}

	@Override
	public String getUniProtId() {
		return uniProtId;
	}

	@Override
	public void setUniProtId(String uniProtId) {
		this.uniProtId = uniProtId;
	}

	@Override
	public String getInsCode() {
		return insCode;
	}

	@Override
	public void setInsCode(String insCode) {
		this.insCode = insCode;
	}

	@Override
	public String getPdbResNum() {
		return pdbResNum;
	}

	@Override
	public void setPdbResNum(String pdbResNum) {
		this.pdbResNum = pdbResNum;
	}

	@Override
	public String toString() {
		StringBuffer  s = new StringBuffer();

		s.append("SeqMisMatchImpl{");
		s.append("seqNum=" );
		s.append(seqNum );
		s.append(", origGroup='" );
		s.append(origGroup + '\'' );
		s.append(", pdbGroup='" );
		s.append(pdbGroup + '\'' );
		s.append(", details='" );
		s.append(details + '\'' );
		s.append(", uniProtId='" );
		s.append(uniProtId + '\'' );
		s.append(", pdbResNum='" );
		s.append(pdbResNum + '\'' );

		if ( insCode == null)
			s.append(", insCode=null ") ;
		else
			s.append(", insCode='" + insCode + '\'') ;

		s.append('}');
		return s.toString();

	}
}
