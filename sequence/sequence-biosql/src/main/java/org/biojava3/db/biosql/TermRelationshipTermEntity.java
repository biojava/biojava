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
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToOne;
import javax.persistence.Table;

/**
 * Entity for TermRelationshipTerm
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "TERM_RELATIONSHIP_TERM")
@NamedQueries({
    @NamedQuery(name = "TermRelationshipTermEntity.findByTermRelationshipId", query = "SELECT t FROM TermRelationshipTermEntity t WHERE t.termRelationshipId = :termRelationshipId")})
public class TermRelationshipTermEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "TERM_RELATIONSHIP_ID", nullable = false)
    private Integer termRelationshipId;
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID")
    @OneToOne
    private TermEntity termId;
    @JoinColumn(name = "TERM_RELATIONSHIP_ID", referencedColumnName = "TERM_RELATIONSHIP_ID", insertable = false, updatable = false)
    @OneToOne
    private TermRelationshipEntity termRelationship;

    public TermRelationshipTermEntity() {
    }

    public TermRelationshipTermEntity(Integer termRelationshipId) {
        this.termRelationshipId = termRelationshipId;
    }

    public Integer getTermRelationshipId() {
        return termRelationshipId;
    }

    public void setTermRelationshipId(Integer termRelationshipId) {
        this.termRelationshipId = termRelationshipId;
    }

    public TermEntity getTermId() {
        return termId;
    }

    public void setTermId(TermEntity termId) {
        this.termId = termId;
    }

    public TermRelationshipEntity getTermRelationship() {
        return termRelationship;
    }

    public void setTermRelationship(TermRelationshipEntity termRelationship) {
        this.termRelationship = termRelationship;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (termRelationshipId != null ? termRelationshipId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof TermRelationshipTermEntity)) {
            return false;
        }
        TermRelationshipTermEntity other = (TermRelationshipTermEntity) object;
        if ((this.termRelationshipId == null && other.termRelationshipId != null) || (this.termRelationshipId != null && !this.termRelationshipId.equals(other.termRelationshipId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TermRelationshipTermEntity[termRelationshipId=" + termRelationshipId + "]";
    }

}
