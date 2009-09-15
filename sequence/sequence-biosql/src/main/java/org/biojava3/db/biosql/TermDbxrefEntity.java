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
 * Entity for TermDbxrefEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "TERM_DBXREF")
@NamedQueries({@NamedQuery(name = "TermDbxrefEntity.findByTermId", query = "SELECT t FROM TermDbxrefEntity t WHERE t.termDbxrefPK.termId = :termId"), 
@NamedQuery(name = "TermDbxrefEntity.findByDbxrefId", query = "SELECT t FROM TermDbxrefEntity t WHERE t.termDbxrefPK.dbxrefId = :dbxrefId"), 
@NamedQuery(name = "TermDbxrefEntity.findByRank", query = "SELECT t FROM TermDbxrefEntity t WHERE t.rank = :rank")})
public class TermDbxrefEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    protected TermDbxrefEntityPK termDbxrefPK;
    @Column(name = "RANK")
    private Integer rank;
    @JoinColumn(name = "DBXREF_ID", referencedColumnName = "DBXREF_ID", insertable = false, updatable = false)
    @ManyToOne
    private DbxrefEntity dbxref;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;

    public TermDbxrefEntity() {
    }

    public TermDbxrefEntity(TermDbxrefEntityPK termDbxrefPK) {
        this.termDbxrefPK = termDbxrefPK;
    }

    public TermDbxrefEntity(int termId, int dbxrefId) {
        this.termDbxrefPK = new TermDbxrefEntityPK(termId, dbxrefId);
    }

    public TermDbxrefEntityPK getTermDbxrefPK() {
        return termDbxrefPK;
    }

    public void setTermDbxrefPK(TermDbxrefEntityPK termDbxrefPK) {
        this.termDbxrefPK = termDbxrefPK;
    }

    public Integer getRank() {
        return rank;
    }

    public void setRank(Integer rank) {
        this.rank = rank;
    }

    public DbxrefEntity getDbxref() {
        return dbxref;
    }

    public void setDbxref(DbxrefEntity dbxref) {
        this.dbxref = dbxref;
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
        hash += (termDbxrefPK != null ? termDbxrefPK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof TermDbxrefEntity)) {
            return false;
        }
        TermDbxrefEntity other = (TermDbxrefEntity) object;
        if ((this.termDbxrefPK == null && other.termDbxrefPK != null) || (this.termDbxrefPK != null && !this.termDbxrefPK.equals(other.termDbxrefPK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TermDbxrefEntity[termDbxrefPK=" + termDbxrefPK + "]";
    }

}
