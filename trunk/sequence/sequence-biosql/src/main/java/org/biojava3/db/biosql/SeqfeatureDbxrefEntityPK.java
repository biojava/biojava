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
 * Embeddable PK for SeqfeatureDbxrefEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class SeqfeatureDbxrefEntityPK implements Serializable {
    private static final long serialVersionUID = -1622478890369860302L;
    @Column(name = "SEQFEATURE_ID", nullable = false)
    private int seqfeatureId;
    @Column(name = "DBXREF_ID", nullable = false)
    private int dbxrefId;

    public SeqfeatureDbxrefEntityPK() {
    }

    public SeqfeatureDbxrefEntityPK(int seqfeatureId, int dbxrefId) {
        this.seqfeatureId = seqfeatureId;
        this.dbxrefId = dbxrefId;
    }

    public int getSeqfeatureId() {
        return seqfeatureId;
    }

    public void setSeqfeatureId(int seqfeatureId) {
        this.seqfeatureId = seqfeatureId;
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
        hash += (int) seqfeatureId;
        hash += (int) dbxrefId;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqfeatureDbxrefEntityPK)) {
            return false;
        }
        SeqfeatureDbxrefEntityPK other = (SeqfeatureDbxrefEntityPK) object;
        if (this.seqfeatureId != other.seqfeatureId) {
            return false;
        }
        if (this.dbxrefId != other.dbxrefId) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.SeqfeatureDbxrefEntityPK[seqfeatureId=" + seqfeatureId + ", dbxrefId=" + dbxrefId + "]";
    }

}
