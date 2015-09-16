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
