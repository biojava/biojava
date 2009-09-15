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
package org.biojava3.db.biosql;

import javax.persistence.Embeddable;

/**
 * Embeddable unique key for SeqFeaturePathEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class SeqFeaturePathEntityUK implements java.io.Serializable {
    private static final long serialVersionUID = -650215803858695626L;
    private int objectSeqfeatureId;
    private int subjectSeqfeatureId;
    private int termId;
    private int distance;
    
    public SeqFeaturePathEntityUK(){}
    
    public SeqFeaturePathEntityUK(int objectSeqfeatureId, int subjectSeqfeatureId, int termId, int distance){
        this.distance = distance;
        this.objectSeqfeatureId = objectSeqfeatureId;
        this.subjectSeqfeatureId = subjectSeqfeatureId;
        this.termId = termId;
    }
    
    @Override
    public int hashCode() {
        int hash = 0;
        hash += (int) getObjectSeqfeatureId();
        hash += (int) getSubjectSeqfeatureId();
        hash += (int) getTermId();
        hash += (int) getDistance();
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqFeaturePathEntityUK)) {
            return false;
        }
        SeqFeaturePathEntityUK other = (SeqFeaturePathEntityUK) object;
        if (this.getSubjectSeqfeatureId() != other.getSubjectSeqfeatureId()) {
            return false;
        }
        if (this.getObjectSeqfeatureId() != other.getObjectSeqfeatureId()) {
            return false;
        }
        if (this.getTermId() != other.getTermId()) {
            return false;
        }
        if (this.getDistance() != other.getDistance()) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.SeqFeaturePathEntityUK[objectSeqfeatureId=" + getObjectSeqfeatureId() +
                "subjectSeqfeatureId=" + getSubjectSeqfeatureId() +
                ", termId=" + getTermId() + ", distance=" + getDistance() + "]";
    }

    public int getObjectSeqfeatureId() {
        return objectSeqfeatureId;
    }

    public void setObjectSeqfeatureId(int objectSeqfeatureId) {
        this.objectSeqfeatureId = objectSeqfeatureId;
    }

    public int getSubjectSeqfeatureId() {
        return subjectSeqfeatureId;
    }

    public void setSubjectSeqfeatureId(int subjectSeqfeatureId) {
        this.subjectSeqfeatureId = subjectSeqfeatureId;
    }

    public int getTermId() {
        return termId;
    }

    public void setTermId(int termId) {
        this.termId = termId;
    }

    public int getDistance() {
        return distance;
    }

    public void setDistance(int distance) {
        this.distance = distance;
    }
}
