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
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToMany;
import javax.persistence.OneToOne;
import javax.persistence.Table;

/**
 * Entity for ReferenceEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "REFERENCE")
@NamedQueries({
    @NamedQuery(name = "ReferenceEntity.findByReferenceId", query = "SELECT r FROM ReferenceEntity r WHERE r.referenceId = :referenceId"), 
    @NamedQuery(name = "ReferenceEntity.findByLocation", query = "SELECT r FROM ReferenceEntity r WHERE r.location = :location"), 
    @NamedQuery(name = "ReferenceEntity.findByTitle", query = "SELECT r FROM ReferenceEntity r WHERE r.title = :title"), 
    @NamedQuery(name = "ReferenceEntity.findByAuthors", query = "SELECT r FROM ReferenceEntity r WHERE r.authors = :authors"), 
    @NamedQuery(name = "ReferenceEntity.findByCrc", query = "SELECT r FROM ReferenceEntity r WHERE r.crc = :crc")})
public class ReferenceEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "REFERENCE_ID", nullable = false)
    private Integer referenceId;
    @Column(name = "LOCATION", nullable = false)
    private String location;
    @Column(name = "TITLE")
    private String title;
    @Column(name = "AUTHORS")
    private String authors;
    @Column(name = "CRC", nullable = false)
    private String crc;
    @JoinColumn(name = "DBXREF_ID", referencedColumnName = "DBXREF_ID")
    @OneToOne
    private DbxrefEntity dbxrefId;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "reference")
    private Collection<BioentryReferenceEntity> bioentryReferenceCollection;

    public ReferenceEntity() {
    }

    public ReferenceEntity(Integer referenceId) {
        this.referenceId = referenceId;
    }

    public ReferenceEntity(Integer referenceId, String location, String crc) {
        this.referenceId = referenceId;
        this.location = location;
        this.crc = crc;
    }

    public Integer getReferenceId() {
        return referenceId;
    }

    public void setReferenceId(Integer referenceId) {
        this.referenceId = referenceId;
    }

    public String getLocation() {
        return location;
    }

    public void setLocation(String location) {
        this.location = location;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    public String getAuthors() {
        return authors;
    }

    public void setAuthors(String authors) {
        this.authors = authors;
    }

    public String getCrc() {
        return crc;
    }

    public void setCrc(String crc) {
        this.crc = crc;
    }

    public DbxrefEntity getDbxrefId() {
        return dbxrefId;
    }

    public void setDbxrefId(DbxrefEntity dbxrefId) {
        this.dbxrefId = dbxrefId;
    }

    public Collection<BioentryReferenceEntity> getBioentryReferenceCollection() {
        return bioentryReferenceCollection;
    }

    public void setBioentryReferenceCollection(Collection<BioentryReferenceEntity> bioentryReferenceCollection) {
        this.bioentryReferenceCollection = bioentryReferenceCollection;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (referenceId != null ? referenceId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof ReferenceEntity)) {
            return false;
        }
        ReferenceEntity other = (ReferenceEntity) object;
        if ((this.referenceId == null && other.referenceId != null) || (this.referenceId != null && !this.referenceId.equals(other.referenceId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.ReferenceEntity[referenceId=" + referenceId + "]";
    }

}
