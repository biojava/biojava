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
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToMany;
import javax.persistence.OneToOne;
import javax.persistence.Table;

/**
 * The persitent entity for a DbxrefEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "DBXREF")
@NamedQueries({
    @NamedQuery(name = "DbxrefEntity.findByDbxrefId", query = "SELECT d FROM DbxrefEntity d WHERE d.dbxrefId = :dbxrefId"), 
    @NamedQuery(name = "DbxrefEntity.findByDbname", query = "SELECT d FROM DbxrefEntity d WHERE d.dbname = :dbname"), 
    @NamedQuery(name = "DbxrefEntity.findByAccession", query = "SELECT d FROM DbxrefEntity d WHERE d.accession = :accession"), 
    @NamedQuery(name = "DbxrefEntity.findByVersion", query = "SELECT d FROM DbxrefEntity d WHERE d.version = :version")})
public class DbxrefEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "DBXREF_ID", nullable = false)
    private Integer dbxrefId;
    @Column(name = "DBNAME", nullable = false)
    private String dbname;
    @Column(name = "ACCESSION", nullable = false)
    private String accession;
    @Column(name = "VERSION", nullable = false)
    private int version;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "dbxref")
    private Collection<TermDbxrefEntity> termDbxrefCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "dbxref")
    private Collection<DbxrefQualifierValueEntity> dbxrefQualifierValueCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "dbxref")
    private Collection<BioentryDbxrefEntity> bioentryDbxrefCollection;
    @OneToMany(mappedBy = "dbxrefId")
    private Collection<LocationEntity> locationCollection;
    @OneToOne(cascade = CascadeType.ALL, mappedBy = "dbxrefId")
    private ReferenceEntity reference;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "dbxref")
    private Collection<SeqfeatureDbxrefEntity> seqfeatureDbxrefCollection;

    public DbxrefEntity() {
    }

    public DbxrefEntity(Integer dbxrefId) {
        this.dbxrefId = dbxrefId;
    }

    public DbxrefEntity(Integer dbxrefId, String dbname, String accession, int version) {
        this.dbxrefId = dbxrefId;
        this.dbname = dbname;
        this.accession = accession;
        this.version = version;
    }

    public Integer getDbxrefId() {
        return dbxrefId;
    }

    public void setDbxrefId(Integer dbxrefId) {
        this.dbxrefId = dbxrefId;
    }

    public String getDbname() {
        return dbname;
    }

    public void setDbname(String dbname) {
        this.dbname = dbname;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public int getVersion() {
        return version;
    }

    public void setVersion(int version) {
        this.version = version;
    }

    public Collection<TermDbxrefEntity> getTermDbxrefCollection() {
        return termDbxrefCollection;
    }

    public void setTermDbxrefCollection(Collection<TermDbxrefEntity> termDbxrefCollection) {
        this.termDbxrefCollection = termDbxrefCollection;
    }

    public Collection<DbxrefQualifierValueEntity> getDbxrefQualifierValueCollection() {
        return dbxrefQualifierValueCollection;
    }

    public void setDbxrefQualifierValueCollection(Collection<DbxrefQualifierValueEntity> dbxrefQualifierValueCollection) {
        this.dbxrefQualifierValueCollection = dbxrefQualifierValueCollection;
    }

    public Collection<BioentryDbxrefEntity> getBioentryDbxrefCollection() {
        return bioentryDbxrefCollection;
    }

    public void setBioentryDbxrefCollection(Collection<BioentryDbxrefEntity> bioentryDbxrefCollection) {
        this.bioentryDbxrefCollection = bioentryDbxrefCollection;
    }

    public Collection<LocationEntity> getLocationCollection() {
        return locationCollection;
    }

    public void setLocationCollection(Collection<LocationEntity> locationCollection) {
        this.locationCollection = locationCollection;
    }

    public ReferenceEntity getReference() {
        return reference;
    }

    public void setReference(ReferenceEntity reference) {
        this.reference = reference;
    }

    public Collection<SeqfeatureDbxrefEntity> getSeqfeatureDbxrefCollection() {
        return seqfeatureDbxrefCollection;
    }

    public void setSeqfeatureDbxrefCollection(Collection<SeqfeatureDbxrefEntity> seqfeatureDbxrefCollection) {
        this.seqfeatureDbxrefCollection = seqfeatureDbxrefCollection;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (dbxrefId != null ? dbxrefId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof DbxrefEntity)) {
            return false;
        }
        DbxrefEntity other = (DbxrefEntity) object;
        if ((this.dbxrefId == null && other.dbxrefId != null) || (this.dbxrefId != null && !this.dbxrefId.equals(other.dbxrefId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.DbxrefEntity[dbxrefId=" + dbxrefId + "]";
    }

}
