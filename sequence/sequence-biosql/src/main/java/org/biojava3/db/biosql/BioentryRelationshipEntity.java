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
 * The entity bean for a BioentryRelationshipEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "BIOENTRY_RELATIONSHIP")
@NamedQueries({
    @NamedQuery(name = "BioentryRelationshipEntity.findByBioentryRelationshipId", query = "SELECT b FROM BioentryRelationshipEntity b WHERE b.bioentryRelationshipId = :bioentryRelationshipId"),
    @NamedQuery(name = "BioentryRelationshipEntity.findByRank", query = "SELECT b FROM BioentryRelationshipEntity b WHERE b.rank = :rank")})
public class BioentryRelationshipEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "BIOENTRY_RELATIONSHIP_ID", nullable = false)
    private Integer bioentryRelationshipId;
    @Column(name = "RANK")
    private Integer rank;
    @JoinColumn(name = "SUBJECT_BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID")
    @ManyToOne
    private BioentryEntity subjectBioentryId;
    @JoinColumn(name = "OBJECT_BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID")
    @ManyToOne
    private BioentryEntity objectBioentryId;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity termId;

    public BioentryRelationshipEntity() {
    }

    public BioentryRelationshipEntity(Integer bioentryRelationshipId) {
        this.bioentryRelationshipId = bioentryRelationshipId;
    }

    public Integer getBioentryRelationshipId() {
        return bioentryRelationshipId;
    }

    public void setBioentryRelationshipId(Integer bioentryRelationshipId) {
        this.bioentryRelationshipId = bioentryRelationshipId;
    }

    public Integer getRank() {
        return rank;
    }

    public void setRank(Integer rank) {
        this.rank = rank;
    }

    public BioentryEntity getSubjectBioentryId() {
        return subjectBioentryId;
    }

    public void setSubjectBioentryId(BioentryEntity subjectBioentryId) {
        this.subjectBioentryId = subjectBioentryId;
    }

    public BioentryEntity getObjectBioentryId() {
        return objectBioentryId;
    }

    public void setObjectBioentryId(BioentryEntity objectBioentryId) {
        this.objectBioentryId = objectBioentryId;
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
        hash += (bioentryRelationshipId != null ? bioentryRelationshipId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryRelationshipEntity)) {
            return false;
        }
        BioentryRelationshipEntity other = (BioentryRelationshipEntity) object;
        if ((this.bioentryRelationshipId == null && other.bioentryRelationshipId != null) || (this.bioentryRelationshipId != null && !this.bioentryRelationshipId.equals(other.bioentryRelationshipId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryRelationshipEntity[bioentryRelationshipId=" + bioentryRelationshipId + "]";
    }

}
