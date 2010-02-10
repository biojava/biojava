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
import javax.persistence.Embeddable;

/**
 * The embedded PK for BioentryReferenceEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class BioentryReferenceEntityPK implements Serializable {
    private static final long serialVersionUID = -1570673432209440298L;
    @Column(name = "BIOENTRY_ID", nullable = false)
    private int bioentryId;
    @Column(name = "REFERENCE_ID", nullable = false)
    private int referenceId;
    @Column(name = "RANK", nullable = false)
    private int rank;

    public BioentryReferenceEntityPK() {
    }

    public BioentryReferenceEntityPK(int bioentryId, int referenceId, int rank) {
        this.bioentryId = bioentryId;
        this.referenceId = referenceId;
        this.rank = rank;
    }

    public int getBioentryId() {
        return bioentryId;
    }

    public void setBioentryId(int bioentryId) {
        this.bioentryId = bioentryId;
    }

    public int getReferenceId() {
        return referenceId;
    }

    public void setReferenceId(int referenceId) {
        this.referenceId = referenceId;
    }

    public int getRank() {
        return rank;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (int) bioentryId;
        hash += (int) referenceId;
        hash += (int) rank;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryReferenceEntityPK)) {
            return false;
        }
        BioentryReferenceEntityPK other = (BioentryReferenceEntityPK) object;
        if (this.bioentryId != other.bioentryId) {
            return false;
        }
        if (this.referenceId != other.referenceId) {
            return false;
        }
        if (this.rank != other.rank) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryReferencePK[bioentryId=" + bioentryId + ", referenceId=" + referenceId + ", rank=" + rank + "]";
    }

}
