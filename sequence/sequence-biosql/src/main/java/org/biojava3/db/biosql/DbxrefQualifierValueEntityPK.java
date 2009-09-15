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
 * The embeddable primary key for DbxrefQualifierValueEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class DbxrefQualifierValueEntityPK implements Serializable {
    private static final long serialVersionUID = 896740444685751623L;
    @Column(name = "DBXREF_ID", nullable = false)
    private int dbxrefId;
    @Column(name = "TERM_ID", nullable = false)
    private int termId;
    @Column(name = "RANK", nullable = false)
    private int rank;

    public DbxrefQualifierValueEntityPK() {
    }

    public DbxrefQualifierValueEntityPK(int dbxrefId, int termId, int rank) {
        this.dbxrefId = dbxrefId;
        this.termId = termId;
        this.rank = rank;
    }

    public int getDbxrefId() {
        return dbxrefId;
    }

    public void setDbxrefId(int dbxrefId) {
        this.dbxrefId = dbxrefId;
    }

    public int getTermId() {
        return termId;
    }

    public void setTermId(int termId) {
        this.termId = termId;
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
        hash += (int) dbxrefId;
        hash += (int) termId;
        hash += (int) rank;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof DbxrefQualifierValueEntityPK)) {
            return false;
        }
        DbxrefQualifierValueEntityPK other = (DbxrefQualifierValueEntityPK) object;
        if (this.dbxrefId != other.dbxrefId) {
            return false;
        }
        if (this.termId != other.termId) {
            return false;
        }
        if (this.rank != other.rank) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.DbxrefQualifierValueEntityPK[dbxrefId=" + dbxrefId + ", termId=" + termId + ", rank=" + rank + "]";
    }

}
