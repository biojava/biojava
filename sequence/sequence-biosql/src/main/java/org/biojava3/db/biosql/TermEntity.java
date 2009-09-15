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
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToMany;
import javax.persistence.OneToOne;
import javax.persistence.Table;

/**
 * Entity for TermEntity
 * @author Mark Schreiber
 */
@Entity
@Table(name = "TERM")
@NamedQueries({@NamedQuery(name = "TermEntity.findByTermId", query = "SELECT t FROM TermEntity t WHERE t.termId = :termId"), @NamedQuery(name = "TermEntity.findByName", query = "SELECT t FROM TermEntity t WHERE t.name = :name"), @NamedQuery(name = "TermEntity.findByDefinition", query = "SELECT t FROM TermEntity t WHERE t.definition = :definition"), @NamedQuery(name = "TermEntity.findByIdentifier", query = "SELECT t FROM TermEntity t WHERE t.identifier = :identifier"), @NamedQuery(name = "TermEntity.findByIsObsolete", query = "SELECT t FROM TermEntity t WHERE t.isObsolete = :isObsolete")})
public class TermEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "TERM_ID", nullable = false)
    private Integer termId;
    @Column(name = "NAME", nullable = false)
    private String name;
    @Column(name = "DEFINITION")
    private String definition;
    @Column(name = "IDENTIFIER", nullable = false)
    private String identifier;
    @Column(name = "IS_OBSOLETE", nullable = false)
    private char isObsolete;
    @JoinColumn(name = "ONTOLOGY_ID", referencedColumnName = "ONTOLOGY_ID")
    @ManyToOne
    private OntologyEntity ontologyId;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "termId")
    private Collection<SeqfeatureRelationshipEntity> seqfeatureRelationshipCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "term")
    private Collection<TermDbxrefEntity> termDbxrefCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "objectTermId")
    private Collection<TermPathEntity> termPathCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "predicateTermId")
    private Collection<TermPathEntity> termPathCollection1;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "subjectTermId")
    private Collection<TermPathEntity> termPathCollection2;
    @OneToOne(cascade = CascadeType.ALL, mappedBy = "termId")
    private TermRelationshipTermEntity termRelationshipTerm;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "term")
    private Collection<DbxrefQualifierValueEntity> dbxrefQualifierValueCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "term")
    private Collection<TermSynonymEntity> termSynonymCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "objectTermId")
    private Collection<TermRelationshipEntity> termRelationshipCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "predicateTermId")
    private Collection<TermRelationshipEntity> termRelationshipCollection1;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "subjectTermId")
    private Collection<TermRelationshipEntity> termRelationshipCollection2;
    @OneToMany(mappedBy = "termId")
    private Collection<LocationEntity> locationCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "term")
    private Collection<SeqfeatureQualifierValueEntity> seqfeatureQualifierValueCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "sourceTermId")
    private Collection<SeqfeatureEntity> seqfeatureCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "typeTermId")
    private Collection<SeqfeatureEntity> seqfeatureCollection1;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "termId")
    private Collection<BioentryRelationshipEntity> bioentryRelationshipCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "term")
    private Collection<LocationQualifierValueEntity> locationQualifierValueCollection;

    public TermEntity() {
    }

    public TermEntity(Integer termId) {
        this.termId = termId;
    }

    public TermEntity(Integer termId, String name, String identifier, char isObsolete) {
        this.termId = termId;
        this.name = name;
        this.identifier = identifier;
        this.isObsolete = isObsolete;
    }

    public Integer getTermId() {
        return termId;
    }

    public void setTermId(Integer termId) {
        this.termId = termId;
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

    public String getIdentifier() {
        return identifier;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    public char getIsObsolete() {
        return isObsolete;
    }

    public void setIsObsolete(char isObsolete) {
        this.isObsolete = isObsolete;
    }

    public OntologyEntity getOntologyId() {
        return ontologyId;
    }

    public void setOntologyId(OntologyEntity ontologyId) {
        this.ontologyId = ontologyId;
    }

    public Collection<SeqfeatureRelationshipEntity> getSeqfeatureRelationshipCollection() {
        return seqfeatureRelationshipCollection;
    }

    public void setSeqfeatureRelationshipCollection(Collection<SeqfeatureRelationshipEntity> seqfeatureRelationshipCollection) {
        this.seqfeatureRelationshipCollection = seqfeatureRelationshipCollection;
    }

    public Collection<TermDbxrefEntity> getTermDbxrefCollection() {
        return termDbxrefCollection;
    }

    public void setTermDbxrefCollection(Collection<TermDbxrefEntity> termDbxrefCollection) {
        this.termDbxrefCollection = termDbxrefCollection;
    }

    public Collection<TermPathEntity> getTermPathCollection() {
        return termPathCollection;
    }

    public void setTermPathCollection(Collection<TermPathEntity> termPathCollection) {
        this.termPathCollection = termPathCollection;
    }

    public Collection<TermPathEntity> getTermPathCollection1() {
        return termPathCollection1;
    }

    public void setTermPathCollection1(Collection<TermPathEntity> termPathCollection1) {
        this.termPathCollection1 = termPathCollection1;
    }

    public Collection<TermPathEntity> getTermPathCollection2() {
        return termPathCollection2;
    }

    public void setTermPathCollection2(Collection<TermPathEntity> termPathCollection2) {
        this.termPathCollection2 = termPathCollection2;
    }

    public TermRelationshipTermEntity getTermRelationshipTerm() {
        return termRelationshipTerm;
    }

    public void setTermRelationshipTerm(TermRelationshipTermEntity termRelationshipTerm) {
        this.termRelationshipTerm = termRelationshipTerm;
    }

    public Collection<DbxrefQualifierValueEntity> getDbxrefQualifierValueCollection() {
        return dbxrefQualifierValueCollection;
    }

    public void setDbxrefQualifierValueCollection(Collection<DbxrefQualifierValueEntity> dbxrefQualifierValueCollection) {
        this.dbxrefQualifierValueCollection = dbxrefQualifierValueCollection;
    }

    public Collection<TermSynonymEntity> getTermSynonymCollection() {
        return termSynonymCollection;
    }

    public void setTermSynonymCollection(Collection<TermSynonymEntity> termSynonymCollection) {
        this.termSynonymCollection = termSynonymCollection;
    }

    public Collection<TermRelationshipEntity> getTermRelationshipCollection() {
        return termRelationshipCollection;
    }

    public void setTermRelationshipCollection(Collection<TermRelationshipEntity> termRelationshipCollection) {
        this.termRelationshipCollection = termRelationshipCollection;
    }

    public Collection<TermRelationshipEntity> getTermRelationshipCollection1() {
        return termRelationshipCollection1;
    }

    public void setTermRelationshipCollection1(Collection<TermRelationshipEntity> termRelationshipCollection1) {
        this.termRelationshipCollection1 = termRelationshipCollection1;
    }

    public Collection<TermRelationshipEntity> getTermRelationshipCollection2() {
        return termRelationshipCollection2;
    }

    public void setTermRelationshipCollection2(Collection<TermRelationshipEntity> termRelationshipCollection2) {
        this.termRelationshipCollection2 = termRelationshipCollection2;
    }

    public Collection<LocationEntity> getLocationCollection() {
        return locationCollection;
    }

    public void setLocationCollection(Collection<LocationEntity> locationCollection) {
        this.locationCollection = locationCollection;
    }

    public Collection<SeqfeatureQualifierValueEntity> getSeqfeatureQualifierValueCollection() {
        return seqfeatureQualifierValueCollection;
    }

    public void setSeqfeatureQualifierValueCollection(Collection<SeqfeatureQualifierValueEntity> seqfeatureQualifierValueCollection) {
        this.seqfeatureQualifierValueCollection = seqfeatureQualifierValueCollection;
    }

    public Collection<SeqfeatureEntity> getSeqfeatureCollection() {
        return seqfeatureCollection;
    }

    public void setSeqfeatureCollection(Collection<SeqfeatureEntity> seqfeatureCollection) {
        this.seqfeatureCollection = seqfeatureCollection;
    }

    public Collection<SeqfeatureEntity> getSeqfeatureCollection1() {
        return seqfeatureCollection1;
    }

    public void setSeqfeatureCollection1(Collection<SeqfeatureEntity> seqfeatureCollection1) {
        this.seqfeatureCollection1 = seqfeatureCollection1;
    }

    public Collection<BioentryRelationshipEntity> getBioentryRelationshipCollection() {
        return bioentryRelationshipCollection;
    }

    public void setBioentryRelationshipCollection(Collection<BioentryRelationshipEntity> bioentryRelationshipCollection) {
        this.bioentryRelationshipCollection = bioentryRelationshipCollection;
    }

    public Collection<LocationQualifierValueEntity> getLocationQualifierValueCollection() {
        return locationQualifierValueCollection;
    }

    public void setLocationQualifierValueCollection(Collection<LocationQualifierValueEntity> locationQualifierValueCollection) {
        this.locationQualifierValueCollection = locationQualifierValueCollection;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (termId != null ? termId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof TermEntity)) {
            return false;
        }
        TermEntity other = (TermEntity) object;
        if ((this.termId == null && other.termId != null) || (this.termId != null && !this.termId.equals(other.termId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TermEntity[termId=" + termId + "]";
    }

}
