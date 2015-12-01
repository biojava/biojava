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

    public Integer getSeqNum() ;

    public void setSeqNum(Integer seqNum) ;

    public String getOrigGroup() ;

    public void setOrigGroup(String origGroup);

    public String getPdbGroup() ;

    public void setPdbGroup(String pdbGroup) ;

    public String getDetails() ;

    public void setDetails(String details);
    public String getUniProtId() ;

    public void setUniProtId(String uniProtId) ;

    public String getInsCode() ;

    public void setInsCode(String insCode) ;

    public String getPdbResNum() ;

    public void setPdbResNum(String pdbResNum) ;

    public String toString();
}
