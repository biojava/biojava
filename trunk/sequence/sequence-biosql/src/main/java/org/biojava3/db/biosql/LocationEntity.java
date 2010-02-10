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
 * Entity bean for LocationEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "LOCATION")
@NamedQueries({
    @NamedQuery(name = "LocationEntity.findByLocationId", query = "SELECT l FROM LocationEntity l WHERE l.locationId = :locationId"), 
    @NamedQuery(name = "LocationEntity.findByStartPos", query = "SELECT l FROM LocationEntity l WHERE l.startPos = :startPos"), 
    @NamedQuery(name = "LocationEntity.findByEndPos", query = "SELECT l FROM LocationEntity l WHERE l.endPos = :endPos"), 
    @NamedQuery(name = "LocationEntity.findByStrand", query = "SELECT l FROM LocationEntity l WHERE l.strand = :strand"), 
    @NamedQuery(name = "LocationEntity.findByRank", query = "SELECT l FROM LocationEntity l WHERE l.rank = :rank")})
public class LocationEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "LOCATION_ID", nullable = false)
    private Integer locationId;
    @Column(name = "START_POS")
    private Integer startPos;
    @Column(name = "END_POS")
    private Integer endPos;
    @Column(name = "STRAND", nullable = false)
    private int strand;
    @Column(name = "RANK", nullable = false)
    private int rank;
    @JoinColumn(name = "DBXREF_ID", referencedColumnName = "DBXREF_ID")
    @ManyToOne
    private DbxrefEntity dbxrefId;
    @JoinColumn(name = "SEQFEATURE_ID", referencedColumnName = "SEQFEATURE_ID")
    @ManyToOne
    private SeqfeatureEntity seqfeatureId;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity termId;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "location")
    private Collection<LocationQualifierValueEntity> locationQualifierValueCollection;

    public LocationEntity() {
    }

    public LocationEntity(Integer locationId) {
        this.locationId = locationId;
    }

    public LocationEntity(Integer locationId, int strand, int rank) {
        this.locationId = locationId;
        this.strand = strand;
        this.rank = rank;
    }

    public Integer getLocationId() {
        return locationId;
    }

    public void setLocationId(Integer locationId) {
        this.locationId = locationId;
    }

    public Integer getStartPos() {
        return startPos;
    }

    public void setStartPos(Integer startPos) {
        this.startPos = startPos;
    }

    public Integer getEndPos() {
        return endPos;
    }

    public void setEndPos(Integer endPos) {
        this.endPos = endPos;
    }

    public int getStrand() {
        return strand;
    }

    public void setStrand(int strand) {
        this.strand = strand;
    }

    public int getRank() {
        return rank;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public DbxrefEntity getDbxrefId() {
        return dbxrefId;
    }

    public void setDbxrefId(DbxrefEntity dbxrefId) {
        this.dbxrefId = dbxrefId;
    }

    public SeqfeatureEntity getSeqfeatureId() {
        return seqfeatureId;
    }

    public void setSeqfeatureId(SeqfeatureEntity seqfeatureId) {
        this.seqfeatureId = seqfeatureId;
    }

    public TermEntity getTermId() {
        return termId;
    }

    public void setTermId(TermEntity termId) {
        this.termId = termId;
    }

    public Collection<LocationQualifierValueEntity> getLocationQualifierValueCollection() {
        return locationQualifierValueCollection;
    }

    public void setLocationQualifierValueCollection(Collection<LocationQualifierValueEntity> locationQualifierValueCollection) {
        this.locationQualifierValueCollection = locationQualifierValueCollection;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (locationId != null ? locationId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof LocationEntity)) {
            return false;
        }
        LocationEntity other = (LocationEntity) object;
        if ((this.locationId == null && other.locationId != null) || (this.locationId != null && !this.locationId.equals(other.locationId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.LocationEntity[locationId=" + locationId + "]";
    }

}
