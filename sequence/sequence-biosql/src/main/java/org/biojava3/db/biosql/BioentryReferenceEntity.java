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
 * The entity bean that backs a BioentryReference implementation
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "BIOENTRY_REFERENCE")
@NamedQueries({ 
    @NamedQuery(name = "BioentryReferenceEntity.findByStartPos", query = "SELECT b FROM BioentryReferenceEntity b WHERE b.startPos = :startPos"), 
    @NamedQuery(name = "BioentryReferenceEntity.findByEndPos", query = "SELECT b FROM BioentryReferenceEntity b WHERE b.endPos = :endPos"), 
    @NamedQuery(name = "BioentryReferenceEntity.findByRank", query = "SELECT b FROM BioentryReferenceEntity b WHERE b.bioentryReferencePK.rank = :rank")})
public class BioentryReferenceEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected BioentryReferenceEntityPK bioentryReferencePK;
    @Column(name = "START_POS")
    private Integer startPos;
    @Column(name = "END_POS")
    private Integer endPos;
    @JoinColumn(name = "BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID", insertable = false, updatable = false)
    @ManyToOne
    private BioentryEntity bioentry;
    @JoinColumn(name = "REFERENCE_ID", referencedColumnName = "REFERENCE_ID", insertable = false, updatable = false)
    @ManyToOne
    private ReferenceEntity reference;

    public BioentryReferenceEntity() {
    }

    public BioentryReferenceEntity(BioentryReferenceEntityPK bioentryReferencePK) {
        this.bioentryReferencePK = bioentryReferencePK;
    }

    public BioentryReferenceEntity(int bioentryId, int referenceId, int rank) {
        this.bioentryReferencePK = new BioentryReferenceEntityPK(bioentryId, referenceId, rank);
    }

    public BioentryReferenceEntityPK getBioentryReferencePK() {
        return bioentryReferencePK;
    }

    public void setBioentryReferencePK(BioentryReferenceEntityPK bioentryReferencePK) {
        this.bioentryReferencePK = bioentryReferencePK;
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

    public BioentryEntity getBioentry() {
        return bioentry;
    }

    public void setBioentry(BioentryEntity bioentry) {
        this.bioentry = bioentry;
    }

    public ReferenceEntity getReference() {
        return reference;
    }

    public void setReference(ReferenceEntity reference) {
        this.reference = reference;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (bioentryReferencePK != null ? bioentryReferencePK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryReferenceEntity)) {
            return false;
        }
        BioentryReferenceEntity other = (BioentryReferenceEntity) object;
        if ((this.bioentryReferencePK == null && other.bioentryReferencePK != null) || (this.bioentryReferencePK != null && !this.bioentryReferencePK.equals(other.bioentryReferencePK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryReferenceEntity[bioentryReferencePK=" + bioentryReferencePK + "]";
    }

}
