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
 * Entity for SeqfeatureDbxref
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "SEQFEATURE_DBXREF")
@NamedQueries({
    @NamedQuery(name = "SeqfeatureDbxrefEntity.findBySeqfeatureId", query = "SELECT s FROM SeqfeatureDbxrefEntity s WHERE s.seqfeatureDbxrefPK.seqfeatureId = :seqfeatureId"), 
    @NamedQuery(name = "SeqfeatureDbxrefEntity.findByDbxrefId", query = "SELECT s FROM SeqfeatureDbxrefEntity s WHERE s.seqfeatureDbxrefPK.dbxrefId = :dbxrefId"), 
    @NamedQuery(name = "SeqfeatureDbxrefEntity.findByRank", query = "SELECT s FROM SeqfeatureDbxrefEntity s WHERE s.rank = :rank")})
public class SeqfeatureDbxrefEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected SeqfeatureDbxrefEntityPK seqfeatureDbxrefPK;
    @Column(name = "RANK")
    private Integer rank;
    @JoinColumn(name = "DBXREF_ID", referencedColumnName = "DBXREF_ID", insertable = false, updatable = false)
    @ManyToOne
    private DbxrefEntity dbxref;
    @JoinColumn(name = "SEQFEATURE_ID", referencedColumnName = "SEQFEATURE_ID", insertable = false, updatable = false)
    @ManyToOne
    private SeqfeatureEntity seqfeature;

    public SeqfeatureDbxrefEntity() {
    }

    public SeqfeatureDbxrefEntity(SeqfeatureDbxrefEntityPK seqfeatureDbxrefPK) {
        this.seqfeatureDbxrefPK = seqfeatureDbxrefPK;
    }

    public SeqfeatureDbxrefEntity(int seqfeatureId, int dbxrefId) {
        this.seqfeatureDbxrefPK = new SeqfeatureDbxrefEntityPK(seqfeatureId, dbxrefId);
    }

    public SeqfeatureDbxrefEntityPK getSeqfeatureDbxrefPK() {
        return seqfeatureDbxrefPK;
    }

    public void setSeqfeatureDbxrefPK(SeqfeatureDbxrefEntityPK seqfeatureDbxrefPK) {
        this.seqfeatureDbxrefPK = seqfeatureDbxrefPK;
    }

    public Integer getRank() {
        return rank;
    }

    public void setRank(Integer rank) {
        this.rank = rank;
    }

    public DbxrefEntity getDbxref() {
        return dbxref;
    }

    public void setDbxref(DbxrefEntity dbxref) {
        this.dbxref = dbxref;
    }

    public SeqfeatureEntity getSeqfeature() {
        return seqfeature;
    }

    public void setSeqfeature(SeqfeatureEntity seqfeature) {
        this.seqfeature = seqfeature;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (seqfeatureDbxrefPK != null ? seqfeatureDbxrefPK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqfeatureDbxrefEntity)) {
            return false;
        }
        SeqfeatureDbxrefEntity other = (SeqfeatureDbxrefEntity) object;
        if ((this.seqfeatureDbxrefPK == null && other.seqfeatureDbxrefPK != null) || (this.seqfeatureDbxrefPK != null && !this.seqfeatureDbxrefPK.equals(other.seqfeatureDbxrefPK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.SeqfeatureDbxrefEntity[seqfeatureDbxrefPK=" + seqfeatureDbxrefPK + "]";
    }

}
