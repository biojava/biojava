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

/**
 * The entity bean that backs a BioentryQualifierValueEntity implementation
 * @author Mark Schreiber
 * @since v3
 */
@Entity
public class BioentryQualifierValueEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    private BioentryQualifierValueEntityUK id;
    @JoinColumn(name = "BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID", insertable = false, updatable = false)
    @ManyToOne
    private BioentryEntity bioentry;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;
    @Column(name = "RANK")
    private Integer rank;
    
    @Column(name = "VALUE")
    private String value;
    

    public BioentryQualifierValueEntity() {
    }

    public BioentryQualifierValueEntity(BioentryQualifierValueEntityUK id) {
        this.id = id;
    }

    public BioentryQualifierValueEntity(int bioentryId, int termId, int rank) {
        this.id = new BioentryQualifierValueEntityUK(bioentryId, termId, rank);
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (getId() != null ? getId().hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryQualifierValueEntity)) {
            return false;
        }
        BioentryQualifierValueEntity other = (BioentryQualifierValueEntity) object;
        if ((this.getId() == null && other.getId() != null) || (this.getId() != null && !this.id.equals(other.id))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryQualifierValue[id=" + getId() + "]";
    }

    public BioentryQualifierValueEntityUK getId() {
        return id;
    }

    public void setId(BioentryQualifierValueEntityUK id) {
        this.id = id;
    }

    public BioentryEntity getBioentry() {
        return bioentry;
    }

    public void setBioentry(BioentryEntity bioentry) {
        this.bioentry = bioentry;
    }

    public TermEntity getTerm() {
        return term;
    }

    public void setTerm(TermEntity term) {
        this.term = term;
    }

    public Integer getRank() {
        return rank;
    }

    public void setRank(Integer rank) {
        this.rank = rank;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }

}
