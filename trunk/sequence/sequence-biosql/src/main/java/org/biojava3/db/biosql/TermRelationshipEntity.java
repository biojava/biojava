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
import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToOne;
import javax.persistence.Table;

/**
 * Entity for TermRelationshipEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "TERM_RELATIONSHIP")
@NamedQueries({
    @NamedQuery(name = "TermRelationshipEntity.findByTermRelationshipId", query = "SELECT t FROM TermRelationshipEntity t WHERE t.termRelationshipId = :termRelationshipId")})
public class TermRelationshipEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "TERM_RELATIONSHIP_ID", nullable = false)
    private Integer termRelationshipId;
    @OneToOne(cascade = CascadeType.ALL, mappedBy = "termRelationship")
    private TermRelationshipTermEntity termRelationshipTerm;
    @JoinColumn(name = "ONTOLOGY_ID", referencedColumnName = "ONTOLOGY_ID")
    @ManyToOne
    private OntologyEntity ontologyId;
    @JoinColumn(name = "OBJECT_TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity objectTermId;
    @JoinColumn(name = "PREDICATE_TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity predicateTermId;
    @JoinColumn(name = "SUBJECT_TERM_ID", referencedColumnName = "TERM_ID")
    @ManyToOne
    private TermEntity subjectTermId;

    public TermRelationshipEntity() {
    }

    public TermRelationshipEntity(Integer termRelationshipId) {
        this.termRelationshipId = termRelationshipId;
    }

    public Integer getTermRelationshipId() {
        return termRelationshipId;
    }

    public void setTermRelationshipId(Integer termRelationshipId) {
        this.termRelationshipId = termRelationshipId;
    }

    public TermRelationshipTermEntity getTermRelationshipTerm() {
        return termRelationshipTerm;
    }

    public void setTermRelationshipTerm(TermRelationshipTermEntity termRelationshipTerm) {
        this.termRelationshipTerm = termRelationshipTerm;
    }

    public OntologyEntity getOntologyId() {
        return ontologyId;
    }

    public void setOntologyId(OntologyEntity ontologyId) {
        this.ontologyId = ontologyId;
    }

    public TermEntity getObjectTermId() {
        return objectTermId;
    }

    public void setObjectTermId(TermEntity objectTermId) {
        this.objectTermId = objectTermId;
    }

    public TermEntity getPredicateTermId() {
        return predicateTermId;
    }

    public void setPredicateTermId(TermEntity predicateTermId) {
        this.predicateTermId = predicateTermId;
    }

    public TermEntity getSubjectTermId() {
        return subjectTermId;
    }

    public void setSubjectTermId(TermEntity subjectTermId) {
        this.subjectTermId = subjectTermId;
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
        if (!(object instanceof TermRelationshipEntity)) {
            return false;
        }
        TermRelationshipEntity other = (TermRelationshipEntity) object;
        if ((this.termRelationshipId == null && other.termRelationshipId != null) || (this.termRelationshipId != null && !this.termRelationshipId.equals(other.termRelationshipId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TermRelationshipEntity[termRelationshipId=" + termRelationshipId + "]";
    }

}
