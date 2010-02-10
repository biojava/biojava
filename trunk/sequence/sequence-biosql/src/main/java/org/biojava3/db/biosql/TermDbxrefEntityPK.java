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
 * Embeddable PK for TermDbxrefEntityPK
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class TermDbxrefEntityPK implements Serializable {
    private static final long serialVersionUID = 6529798870152902086L;
    @Column(name = "TERM_ID", nullable = false)
    private int termId;
    @Column(name = "DBXREF_ID", nullable = false)
    private int dbxrefId;

    public TermDbxrefEntityPK() {
    }

    public TermDbxrefEntityPK(int termId, int dbxrefId) {
        this.termId = termId;
        this.dbxrefId = dbxrefId;
    }

    public int getTermId() {
        return termId;
    }

    public void setTermId(int termId) {
        this.termId = termId;
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
        hash += (int) termId;
        hash += (int) dbxrefId;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof TermDbxrefEntityPK)) {
            return false;
        }
        TermDbxrefEntityPK other = (TermDbxrefEntityPK) object;
        if (this.termId != other.termId) {
            return false;
        }
        if (this.dbxrefId != other.dbxrefId) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TermDbxrefEntityPK[termId=" + termId + ", dbxrefId=" + dbxrefId + "]";
    }

}
