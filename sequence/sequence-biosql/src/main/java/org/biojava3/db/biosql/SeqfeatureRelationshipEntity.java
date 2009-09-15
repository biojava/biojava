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
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.Table;

/**
 * Entity for SeqfeatureRelationshipEntity
 * @author Mark Schreiber
 */
@Entity
@Table(name = "SEQFEATURE_RELATIONSHIP")
@NamedQueries({
    @NamedQuery(name = "SeqfeatureRelationshipEntity.findBySeqfeatureRelationshipId", query = "SELECT s FROM SeqfeatureRelationshipEntity s WHERE s.seqfeatureRelationshipId = :seqfeatureRelationshipId"), 
    @NamedQuery(name = "SeqfeatureRelationshipEntity.findByRank", query = "SELECT s FROM SeqfeatureRelationshipEntity s WHERE s.rank = :rank")})
public class SeqfeatureRelationshipEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "SEQFEATURE_RELATIONSHIP_ID", nullable = false)
    private Integer seqfeatureRelationshipId;
    @Column(name = "RANK")
    private Integer rank;
    @JoinColumn(name = "SUBJECT_SEQFEATURE_ID", referencedColumnName = "SEQFEATURE_ID")
    @ManyToOne
    private SeqfeatureEntity subjectSeqfeatureId;
    @JoinColumn(name = "OBJECT_SEQFEATURE_ID", referencedColumnName = "SEQFEATURE_ID")
    @ManyToOne
    private SeqfeatureEntity objectSeqfeatureId;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity termId;

    public SeqfeatureRelationshipEntity() {
    }

    public SeqfeatureRelationshipEntity(Integer seqfeatureRelationshipId) {
        this.seqfeatureRelationshipId = seqfeatureRelationshipId;
    }

    public Integer getSeqfeatureRelationshipId() {
        return seqfeatureRelationshipId;
    }

    public void setSeqfeatureRelationshipId(Integer seqfeatureRelationshipId) {
        this.seqfeatureRelationshipId = seqfeatureRelationshipId;
    }

    public Integer getRank() {
        return rank;
    }

    public void setRank(Integer rank) {
        this.rank = rank;
    }

    public SeqfeatureEntity getSubjectSeqfeatureId() {
        return subjectSeqfeatureId;
    }

    public void setSubjectSeqfeatureId(SeqfeatureEntity subjectSeqfeatureId) {
        this.subjectSeqfeatureId = subjectSeqfeatureId;
    }

    public SeqfeatureEntity getObjectSeqfeatureId() {
        return objectSeqfeatureId;
    }

    public void setObjectSeqfeatureId(SeqfeatureEntity objectSeqfeatureId) {
        this.objectSeqfeatureId = objectSeqfeatureId;
    }

    public TermEntity getTermId() {
        return termId;
    }

    public void setTermId(TermEntity termId) {
        this.termId = termId;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (seqfeatureRelationshipId != null ? seqfeatureRelationshipId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqfeatureRelationshipEntity)) {
            return false;
        }
        SeqfeatureRelationshipEntity other = (SeqfeatureRelationshipEntity) object;
        if ((this.seqfeatureRelationshipId == null && other.seqfeatureRelationshipId != null) || (this.seqfeatureRelationshipId != null && !this.seqfeatureRelationshipId.equals(other.seqfeatureRelationshipId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.SeqfeatureRelationshipEntity[seqfeatureRelationshipId=" + seqfeatureRelationshipId + "]";
    }

}
