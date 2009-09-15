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
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.Table;

/**
 * The persitent entity bean that backs a BioentryDbxref implementation.
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "BIOENTRY_DBXREF")
@NamedQueries({
    @NamedQuery(name = "BioentryDbxrefEntity.findByBioentryId", query = "SELECT b FROM BioentryDbxrefEntity b WHERE b.bioentryDbxrefPK.bioentryId = :bioentryId"), 
    @NamedQuery(name = "BioentryDbxrefEntity.findByDbxrefId", query = "SELECT b FROM BioentryDbxrefEntity b WHERE b.bioentryDbxrefPK.dbxrefId = :dbxrefId"), 
    @NamedQuery(name = "BioentryDbxrefEntity.findByRank", query = "SELECT b FROM BioentryDbxrefEntity b WHERE b.rank = :rank")})
public class BioentryDbxrefEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected BioentryDbxrefEntityPK bioentryDbxrefPK;
    @Column(name = "RANK")
    private Integer rank;
    @JoinColumn(name = "BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID", insertable = false, updatable = false)
    @ManyToOne
    private BioentryEntity bioentry;
    @JoinColumn(name = "DBXREF_ID", referencedColumnName = "DBXREF_ID", insertable = false, updatable = false)
    @ManyToOne
    private DbxrefEntity dbxref;

    public BioentryDbxrefEntity() {
    }

    public BioentryDbxrefEntity(BioentryDbxrefEntityPK bioentryDbxrefPK) {
        this.bioentryDbxrefPK = bioentryDbxrefPK;
    }

    public BioentryDbxrefEntity(int bioentryId, int dbxrefId) {
        this.bioentryDbxrefPK = new BioentryDbxrefEntityPK(bioentryId, dbxrefId);
    }

    public BioentryDbxrefEntityPK getBioentryDbxrefPK() {
        return bioentryDbxrefPK;
    }

    public void setBioentryDbxrefPK(BioentryDbxrefEntityPK bioentryDbxrefPK) {
        this.bioentryDbxrefPK = bioentryDbxrefPK;
    }

    public Integer getRank() {
        return rank;
    }

    public void setRank(Integer rank) {
        this.rank = rank;
    }

    public BioentryEntity getBioentry() {
        return bioentry;
    }

    public void setBioentry(BioentryEntity bioentry) {
        this.bioentry = bioentry;
    }

    public DbxrefEntity getDbxref() {
        return dbxref;
    }

    public void setDbxref(DbxrefEntity dbxref) {
        this.dbxref = dbxref;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (bioentryDbxrefPK != null ? bioentryDbxrefPK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryDbxrefEntity)) {
            return false;
        }
        BioentryDbxrefEntity other = (BioentryDbxrefEntity) object;
        if ((this.bioentryDbxrefPK == null && other.bioentryDbxrefPK != null) || (this.bioentryDbxrefPK != null && !this.bioentryDbxrefPK.equals(other.bioentryDbxrefPK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryDbxrefEntity[bioentryDbxrefPK=" + bioentryDbxrefPK + "]";
    }

}
