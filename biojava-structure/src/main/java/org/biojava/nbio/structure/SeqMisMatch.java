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

/**
 * Created by andreas on 9/11/15.
 */
public interface SeqMisMatch {

	Integer getSeqNum() ;

	void setSeqNum(Integer seqNum) ;

	String getOrigGroup() ;

	void setOrigGroup(String origGroup);

	String getPdbGroup() ;

	void setPdbGroup(String pdbGroup) ;

	String getDetails() ;

	void setDetails(String details);
	String getUniProtId() ;

	void setUniProtId(String uniProtId) ;

	String getInsCode() ;

	void setInsCode(String insCode) ;

	String getPdbResNum() ;

	void setPdbResNum(String pdbResNum) ;

	@Override
    String toString();
}
