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
 * An embedded PK for BioentryDbXrefEntity. Needed due to the compound key
 * on that table.
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class BioentryDbxrefEntityPK implements Serializable {
    private static final long serialVersionUID = -1577005199688596552L;
    @Column(name = "BIOENTRY_ID", nullable = false)
    private int bioentryId;
    @Column(name = "DBXREF_ID", nullable = false)
    private int dbxrefId;

    public BioentryDbxrefEntityPK() {
    }

    public BioentryDbxrefEntityPK(int bioentryId, int dbxrefId) {
        this.bioentryId = bioentryId;
        this.dbxrefId = dbxrefId;
    }

    public int getBioentryId() {
        return bioentryId;
    }

    public void setBioentryId(int bioentryId) {
        this.bioentryId = bioentryId;
    }

    public int getDbxrefId() {
        return dbxrefId;
    }

    public void setDbxrefId(int dbxrefId) {
        this.dbxrefId = dbxrefId;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (int) bioentryId;
        hash += (int) dbxrefId;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryDbxrefEntityPK)) {
            return false;
        }
        BioentryDbxrefEntityPK other = (BioentryDbxrefEntityPK) object;
        if (this.bioentryId != other.bioentryId) {
            return false;
        }
        if (this.dbxrefId != other.dbxrefId) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryDbxrefPK[bioentryId=" + bioentryId + ", dbxrefId=" + dbxrefId + "]";
    }

}
