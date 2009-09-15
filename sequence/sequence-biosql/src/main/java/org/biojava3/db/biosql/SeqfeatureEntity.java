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
import java.util.Collection;
import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToMany;
import javax.persistence.Table;

/**
 * Entity for SeqfeatureEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "SEQFEATURE")
@NamedQueries({
    @NamedQuery(name = "SeqfeatureEntity.findBySeqfeatureId", query = "SELECT s FROM SeqfeatureEntity s WHERE s.seqfeatureId = :seqfeatureId"), 
    @NamedQuery(name = "SeqfeatureEntity.findByDisplayName", query = "SELECT s FROM SeqfeatureEntity s WHERE s.displayName = :displayName"), 
    @NamedQuery(name = "SeqfeatureEntity.findByRank", query = "SELECT s FROM SeqfeatureEntity s WHERE s.rank = :rank")})
public class SeqfeatureEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "SEQFEATURE_ID", nullable = false)
    private Integer seqfeatureId;
    @Column(name = "DISPLAY_NAME")
    private String displayName;
    @Column(name = "RANK", nullable = false)
    private int rank;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "subjectSeqfeatureId")
    private Collection<SeqfeatureRelationshipEntity> subjectSeqfeatureRelationshipCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "objectSeqfeatureId")
    private Collection<SeqfeatureRelationshipEntity> objectSeqfeatureRelationshipCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "seqfeatureId")
    private Collection<LocationEntity> locationCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "seqfeature")
    private Collection<SeqfeatureDbxrefEntity> seqfeatureDbxrefCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "seqfeature")
    private Collection<SeqfeatureQualifierValueEntity> seqfeatureQualifierValueCollection;
    @JoinColumn(name = "BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID")
    @ManyToOne
    private BioentryEntity bioentryId;
    @JoinColumn(name = "SOURCE_TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity sourceTermId;
    @JoinColumn(name = "TYPE_TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity typeTermId;

    public SeqfeatureEntity() {
    }

    public SeqfeatureEntity(Integer seqfeatureId) {
        this.seqfeatureId = seqfeatureId;
    }

    public SeqfeatureEntity(Integer seqfeatureId, int rank) {
        this.seqfeatureId = seqfeatureId;
        this.rank = rank;
    }

    public Integer getSeqfeatureId() {
        return seqfeatureId;
    }

    public void setSeqfeatureId(Integer seqfeatureId) {
        this.seqfeatureId = seqfeatureId;
    }

    public String getDisplayName() {
        return displayName;
    }

    public void setDisplayName(String displayName) {
        this.displayName = displayName;
    }

    public int getRank() {
        return rank;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public Collection<LocationEntity> getLocationCollection() {
        return locationCollection;
    }

    public void setLocationCollection(Collection<LocationEntity> locationCollection) {
        this.locationCollection = locationCollection;
    }

    public Collection<SeqfeatureDbxrefEntity> getSeqfeatureDbxrefCollection() {
        return seqfeatureDbxrefCollection;
    }

    public void setSeqfeatureDbxrefCollection(Collection<SeqfeatureDbxrefEntity> seqfeatureDbxrefCollection) {
        this.seqfeatureDbxrefCollection = seqfeatureDbxrefCollection;
    }

    public Collection<SeqfeatureQualifierValueEntity> getSeqfeatureQualifierValueCollection() {
        return seqfeatureQualifierValueCollection;
    }

    public void setSeqfeatureQualifierValueCollection(Collection<SeqfeatureQualifierValueEntity> seqfeatureQualifierValueCollection) {
        this.seqfeatureQualifierValueCollection = seqfeatureQualifierValueCollection;
    }

    public BioentryEntity getBioentryId() {
        return bioentryId;
    }

    public void setBioentryId(BioentryEntity bioentryId) {
        this.bioentryId = bioentryId;
    }

    public TermEntity getSourceTermId() {
        return sourceTermId;
    }

    public void setSourceTermId(TermEntity sourceTermId) {
        this.sourceTermId = sourceTermId;
    }

    public TermEntity getTypeTermId() {
        return typeTermId;
    }

    public void setTypeTermId(TermEntity typeTermId) {
        this.typeTermId = typeTermId;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (seqfeatureId != null ? seqfeatureId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqfeatureEntity)) {
            return false;
        }
        SeqfeatureEntity other = (SeqfeatureEntity) object;
        if ((this.seqfeatureId == null && other.seqfeatureId != null) || (this.seqfeatureId != null && !this.seqfeatureId.equals(other.seqfeatureId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.SeqfeatureEntity[seqfeatureId=" + seqfeatureId + "]";
    }

    public Collection<SeqfeatureRelationshipEntity> getSubjectSeqfeatureRelationshipCollection() {
        return subjectSeqfeatureRelationshipCollection;
    }

    public void setSubjectSeqfeatureRelationshipCollection(Collection<SeqfeatureRelationshipEntity> subjectSeqfeatureRelationshipCollection) {
        this.subjectSeqfeatureRelationshipCollection = subjectSeqfeatureRelationshipCollection;
    }

    public Collection<SeqfeatureRelationshipEntity> getObjectSeqfeatureRelationshipCollection() {
        return objectSeqfeatureRelationshipCollection;
    }

    public void setObjectSeqfeatureRelationshipCollection(Collection<SeqfeatureRelationshipEntity> objectSeqfeatureRelationshipCollection) {
        this.objectSeqfeatureRelationshipCollection = objectSeqfeatureRelationshipCollection;
    }

}
