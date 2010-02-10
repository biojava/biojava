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
 * Entity for SeqfeatureQualifierValueEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "SEQFEATURE_QUALIFIER_VALUE")
@NamedQueries({@NamedQuery(name = "SeqfeatureQualifierValueEntity.findBySeqfeatureId", query = "SELECT s FROM SeqfeatureQualifierValueEntity s WHERE s.seqfeatureQualifierValuePK.seqfeatureId = :seqfeatureId"), @NamedQuery(name = "SeqfeatureQualifierValueEntity.findByTermId", query = "SELECT s FROM SeqfeatureQualifierValueEntity s WHERE s.seqfeatureQualifierValuePK.termId = :termId"), @NamedQuery(name = "SeqfeatureQualifierValueEntity.findByRank", query = "SELECT s FROM SeqfeatureQualifierValueEntity s WHERE s.seqfeatureQualifierValuePK.rank = :rank"), @NamedQuery(name = "SeqfeatureQualifierValueEntity.findByValue", query = "SELECT s FROM SeqfeatureQualifierValueEntity s WHERE s.value = :value")})
public class SeqfeatureQualifierValueEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected SeqfeatureQualifierValueEntityPK seqfeatureQualifierValuePK;
    @Column(name = "VALUE", nullable = false)
    private String value;
    @JoinColumn(name = "SEQFEATURE_ID", referencedColumnName = "SEQFEATURE_ID", insertable = false, updatable = false)
    @ManyToOne
    private SeqfeatureEntity seqfeature;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;

    public SeqfeatureQualifierValueEntity() {
    }

    public SeqfeatureQualifierValueEntity(SeqfeatureQualifierValueEntityPK seqfeatureQualifierValuePK) {
        this.seqfeatureQualifierValuePK = seqfeatureQualifierValuePK;
    }

    public SeqfeatureQualifierValueEntity(SeqfeatureQualifierValueEntityPK seqfeatureQualifierValuePK, String value) {
        this.seqfeatureQualifierValuePK = seqfeatureQualifierValuePK;
        this.value = value;
    }

    public SeqfeatureQualifierValueEntity(int seqfeatureId, int termId, int rank) {
        this.seqfeatureQualifierValuePK = new SeqfeatureQualifierValueEntityPK(seqfeatureId, termId, rank);
    }

    public SeqfeatureQualifierValueEntityPK getSeqfeatureQualifierValuePK() {
        return seqfeatureQualifierValuePK;
    }

    public void setSeqfeatureQualifierValuePK(SeqfeatureQualifierValueEntityPK seqfeatureQualifierValuePK) {
        this.seqfeatureQualifierValuePK = seqfeatureQualifierValuePK;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }

    public SeqfeatureEntity getSeqfeature() {
        return seqfeature;
    }

    public void setSeqfeature(SeqfeatureEntity seqfeature) {
        this.seqfeature = seqfeature;
    }

    public TermEntity getTerm() {
        return term;
    }

    public void setTerm(TermEntity term) {
        this.term = term;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (seqfeatureQualifierValuePK != null ? seqfeatureQualifierValuePK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof SeqfeatureQualifierValueEntity)) {
            return false;
        }
        SeqfeatureQualifierValueEntity other = (SeqfeatureQualifierValueEntity) object;
        if ((this.seqfeatureQualifierValuePK == null && other.seqfeatureQualifierValuePK != null) || (this.seqfeatureQualifierValuePK != null && !this.seqfeatureQualifierValuePK.equals(other.seqfeatureQualifierValuePK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.SeqfeatureQualifierValueEntity[seqfeatureQualifierValuePK=" + seqfeatureQualifierValuePK + "]";
    }

}
