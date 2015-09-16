package org.biojava.nbio.structure;

import java.io.Serializable;

/**
 * Created by andreas on 9/11/15.
 */
public class SeqMisMatchImpl implements SeqMisMatch, Serializable{

    Integer seqNum;
    String origGroup;
    String pdbGroup;
    String details;
    String uniProtId;
    String insCode;
    String pdbResNum;

    public Integer getSeqNum() {
        return seqNum;
    }

    public void setSeqNum(Integer seqNum) {
        this.seqNum = seqNum;
    }

    public String getOrigGroup() {
        return origGroup;
    }

    public void setOrigGroup(String origGroup) {
        this.origGroup = origGroup;
    }

    public String getPdbGroup() {
        return pdbGroup;
    }

    public void setPdbGroup(String pdbGroup) {
        this.pdbGroup = pdbGroup;
    }

    public String getDetails() {
        return details;
    }

    public void setDetails(String details) {
        this.details = details;
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public void setUniProtId(String uniProtId) {
        this.uniProtId = uniProtId;
    }

    public String getInsCode() {
        return insCode;
    }

    public void setInsCode(String insCode) {
        this.insCode = insCode;
    }

    public String getPdbResNum() {
        return pdbResNum;
    }

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
