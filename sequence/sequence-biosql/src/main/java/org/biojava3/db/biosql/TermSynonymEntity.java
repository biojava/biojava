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
import javax.persistence.EmbeddedId;
import javax.persistence.Entity;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.Table;

/**
 * Entity for TermSynonym
 * @author Mark
 * @since v3
 */
@Entity
@Table(name = "TERM_SYNONYM")
@NamedQueries({
    @NamedQuery(name = "TermSynonymEntity.findBySynonym", query = "SELECT t FROM TermSynonymEntity t WHERE t.termSynonymPK.synonym = :synonym"), 
    @NamedQuery(name = "TermSynonymEntity.findByTermId", query = "SELECT t FROM TermSynonymEntity t WHERE t.termSynonymPK.termId = :termId")})
public class TermSynonymEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected TermSynonymEntityPK termSynonymPK;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;

    public TermSynonymEntity() {
    }

    public TermSynonymEntity(TermSynonymEntityPK termSynonymPK) {
        this.termSynonymPK = termSynonymPK;
    }

    public TermSynonymEntity(String synonym, int termId) {
        this.termSynonymPK = new TermSynonymEntityPK(synonym, termId);
    }

    public TermSynonymEntityPK getTermSynonymPK() {
        return termSynonymPK;
    }

    public void setTermSynonymPK(TermSynonymEntityPK termSynonymPK) {
        this.termSynonymPK = termSynonymPK;
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
        hash += (termSynonymPK != null ? termSynonymPK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof TermSynonymEntity)) {
            return false;
        }
        TermSynonymEntity other = (TermSynonymEntity) object;
        if ((this.termSynonymPK == null && other.termSynonymPK != null) || (this.termSynonymPK != null && !this.termSynonymPK.equals(other.termSynonymPK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TermSynonymEntity[termSynonymPK=" + termSynonymPK + "]";
    }

}
