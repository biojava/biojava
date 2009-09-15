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

import java.io.Serializable;
import javax.persistence.Column;
import javax.persistence.EmbeddedId;
import javax.persistence.Entity;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;

/**
 * Entity for SeqFeaturePathEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
public class SeqFeaturePathEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    private SeqFeaturePathEntityUK id;
    
    @JoinColumn(name = "OBJECT_SEQFEATURE_ID", referencedColumnName = "OBJECT_SEQFEATURE_ID", insertable = false, updatable = false)
    @ManyToOne
    private SeqfeatureEntity objectSeqfeature;
    @JoinColumn(name = "SUBJECT_SEQFEATURE_ID", referencedColumnName = "SUBJECT_SEQFEATURE_ID", insertable = false, updatable = false)
    @ManyToOne
    private SeqfeatureEntity subjectSeqfeature;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;
    @Column(name="DISTANCE")
    private Integer distance;
    
    public SeqFeaturePathEntity(){}
    public SeqFeaturePathEntity(SeqFeaturePathEntityUK id){
        this.id = id;
    }
    public SeqFeaturePathEntity(int objectSeqfeatureId, int subjectSeqfeatureId, int termId, int distance){
        this.id = new SeqFeaturePathEntityUK(objectSeqfeatureId, subjectSeqfeatureId, termId, distance);
    }


    @Override
    public int hashCode() {
        int hash = 0;
        hash += (id != null ? id.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqFeaturePathEntity)) {
            return false;
        }
        SeqFeaturePathEntity other = (SeqFeaturePathEntity) object;
        if ((this.id == null && other.id != null) || (this.id != null && !this.id.equals(other.id))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.SeqFeaturePathEntity[id=" + id + "]";
    }

    public SeqfeatureEntity getObjectSeqfeature() {
        return objectSeqfeature;
    }

    public void setObjectSeqfeature(SeqfeatureEntity objectSeqfeature) {
        this.objectSeqfeature = objectSeqfeature;
    }

    public SeqfeatureEntity getSubjectSeqfeature() {
        return subjectSeqfeature;
    }

    public void setSubjectSeqfeature(SeqfeatureEntity subjectSeqfeature) {
        this.subjectSeqfeature = subjectSeqfeature;
    }

    public TermEntity getTerm() {
        return term;
    }

    public void setTerm(TermEntity term) {
        this.term = term;
    }

    public Integer getDistance() {
        return distance;
    }

    public void setDistance(Integer distance) {
        this.distance = distance;
    }

}
