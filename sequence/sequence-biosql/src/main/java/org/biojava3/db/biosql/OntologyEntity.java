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
import java.util.Collection;
import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToMany;
import javax.persistence.Table;

/**
 * Entity for OntologyEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "ONTOLOGY")
@NamedQueries({
    @NamedQuery(name = "OntologyEntity.findByOntologyId", query = "SELECT o FROM OntologyEntity o WHERE o.ontologyId = :ontologyId"), 
    @NamedQuery(name = "OntologyEntity.findByName", query = "SELECT o FROM OntologyEntity o WHERE o.name = :name"), 
    @NamedQuery(name = "OntologyEntity.findByDefinition", query = "SELECT o FROM OntologyEntity o WHERE o.definition = :definition")})
public class OntologyEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "ONTOLOGY_ID", nullable = false)
    private Integer ontologyId;
    @Column(name = "NAME", nullable = false)
    private String name;
    @Column(name = "DEFINITION")
    private String definition;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "ontologyId")
    private Collection<TermEntity> termCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "ontologyId")
    private Collection<TermPathEntity> termPathCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "ontologyId")
    private Collection<TermRelationshipEntity> termRelationshipCollection;

    public OntologyEntity() {
    }

    public OntologyEntity(Integer ontologyId) {
        this.ontologyId = ontologyId;
    }

    public OntologyEntity(Integer ontologyId, String name) {
        this.ontologyId = ontologyId;
        this.name = name;
    }

    public Integer getOntologyId() {
        return ontologyId;
    }

    public void setOntologyId(Integer ontologyId) {
        this.ontologyId = ontologyId;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getDefinition() {
        return definition;
    }

    public void setDefinition(String definition) {
        this.definition = definition;
    }

    public Collection<TermEntity> getTermCollection() {
        return termCollection;
    }

    public void setTermCollection(Collection<TermEntity> termCollection) {
        this.termCollection = termCollection;
    }

    public Collection<TermPathEntity> getTermPathCollection() {
        return termPathCollection;
    }

    public void setTermPathCollection(Collection<TermPathEntity> termPathCollection) {
        this.termPathCollection = termPathCollection;
    }

    public Collection<TermRelationshipEntity> getTermRelationshipCollection() {
        return termRelationshipCollection;
    }

    public void setTermRelationshipCollection(Collection<TermRelationshipEntity> termRelationshipCollection) {
        this.termRelationshipCollection = termRelationshipCollection;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (ontologyId != null ? ontologyId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof OntologyEntity)) {
            return false;
        }
        OntologyEntity other = (OntologyEntity) object;
        if ((this.ontologyId == null && other.ontologyId != null) || (this.ontologyId != null && !this.ontologyId.equals(other.ontologyId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.OntologyEntity[ontologyId=" + ontologyId + "]";
    }

}
