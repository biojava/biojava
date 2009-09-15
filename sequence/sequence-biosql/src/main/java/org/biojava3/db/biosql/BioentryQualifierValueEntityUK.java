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

import javax.persistence.Column;
import javax.persistence.Embeddable;

/**
 * The embedded unique key for BioentryQualifierValueEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class BioentryQualifierValueEntityUK implements java.io.Serializable {
    private static final long serialVersionUID = 2581975861292964811L;
    @Column(name = "BIOENTRY_ID", nullable = false)
    private int bioentryId;
    @Column(name = "TERM_ID", nullable = false)
    private int termId;
    @Column(name = "RANK", nullable = false)
    private int rank;
    
    public BioentryQualifierValueEntityUK(){}
    
    public BioentryQualifierValueEntityUK(int bioentryId, int termId, int rank){
        this.bioentryId = bioentryId;
        this.termId = termId;
        this.rank = rank;
    }

    public int getBioentryId() {
        return bioentryId;
    }

    public void setBioentryId(int bioentryId) {
        this.bioentryId = bioentryId;
    }

    public int getRank() {
        return rank;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public int getTermId() {
        return termId;
    }

    public void setTermId(int termId) {
        this.termId = termId;
    }
    @Override
    public int hashCode() {
        int hash = 0;
        hash += (int) bioentryId;
        hash += (int) termId;
        hash += (int) rank;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryQualifierValueEntityUK)) {
            return false;
        }
        BioentryQualifierValueEntityUK other = (BioentryQualifierValueEntityUK) object;
        if (this.bioentryId != other.bioentryId) {
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
        return "org.biosql.entity.BioentryQualifierValueUK[bioentryId=" + bioentryId +
                ", termId=" + termId + ", rank=" + rank + "]";
    }
  
}
