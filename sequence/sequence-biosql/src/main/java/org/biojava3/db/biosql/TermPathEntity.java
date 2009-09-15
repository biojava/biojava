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
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.Table;

/**
 * Entity for TermPathEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "TERM_PATH")
@NamedQueries({
    @NamedQuery(name = "TermPathEntity.findByTermPathId", query = "SELECT t FROM TermPathEntity t WHERE t.termPathId = :termPathId"), 
    @NamedQuery(name = "TermPathEntity.findByDistance", query = "SELECT t FROM TermPathEntity t WHERE t.distance = :distance")})
public class TermPathEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "TERM_PATH_ID", nullable = false)
    private Integer termPathId;
    @Column(name = "DISTANCE", nullable = false)
    private int distance;
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

    public TermPathEntity() {
    }

    public TermPathEntity(Integer termPathId) {
        this.termPathId = termPathId;
    }

    public TermPathEntity(Integer termPathId, int distance) {
        this.termPathId = termPathId;
        this.distance = distance;
    }

    public Integer getTermPathId() {
        return termPathId;
    }

    public void setTermPathId(Integer termPathId) {
        this.termPathId = termPathId;
    }

    public int getDistance() {
        return distance;
    }

    public void setDistance(int distance) {
        this.distance = distance;
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
        hash += (termPathId != null ? termPathId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof TermPathEntity)) {
            return false;
        }
        TermPathEntity other = (TermPathEntity) object;
        if ((this.termPathId == null && other.termPathId != null) || (this.termPathId != null && !this.termPathId.equals(other.termPathId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TermPathEntity[termPathId=" + termPathId + "]";
    }

}
