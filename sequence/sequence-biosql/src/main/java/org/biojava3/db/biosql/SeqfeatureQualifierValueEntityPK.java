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
 * Embeddable PK for SeqfeatureQualifierValueEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class SeqfeatureQualifierValueEntityPK implements Serializable {
    private static final long serialVersionUID = -4873350643917023535L;
    @Column(name = "SEQFEATURE_ID", nullable = false)
    private int seqfeatureId;
    @Column(name = "TERM_ID", nullable = false)
    private int termId;
    @Column(name = "RANK", nullable = false)
    private int rank;

    public SeqfeatureQualifierValueEntityPK() {
    }

    public SeqfeatureQualifierValueEntityPK(int seqfeatureId, int termId, int rank) {
        this.seqfeatureId = seqfeatureId;
        this.termId = termId;
        this.rank = rank;
    }

    public int getSeqfeatureId() {
        return seqfeatureId;
    }

    public void setSeqfeatureId(int seqfeatureId) {
        this.seqfeatureId = seqfeatureId;
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
        hash += (int) seqfeatureId;
        hash += (int) termId;
        hash += (int) rank;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqfeatureQualifierValueEntityPK)) {
            return false;
        }
        SeqfeatureQualifierValueEntityPK other = (SeqfeatureQualifierValueEntityPK) object;
        if (this.seqfeatureId != other.seqfeatureId) {
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
        return "org.biosql.entity.SeqfeatureQualifierValueEntityPK[seqfeatureId=" + seqfeatureId + ", termId=" + termId + ", rank=" + rank + "]";
    }

}
