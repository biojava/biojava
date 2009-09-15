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
import javax.persistence.Table;

/**
 * The entity bean that backs a BioentryPathEntity implementation
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "BIOENTRY_PATH")
public class BioentryPathEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    
    @EmbeddedId
    protected BioentryPathEntityUK id;
    
    @JoinColumn(name = "OBJECT_BIOENTRY_ID", referencedColumnName = "OBJECT_BIOENTRY_ID", insertable = false, updatable = false)
    @ManyToOne
    private BioentryEntity objectBioEntry;
    
    @JoinColumn(name = "SUBJECT_BIOENTRY_ID", referencedColumnName = "SUBJECT_BIOENTRY_ID", insertable = false, updatable = false)
    @ManyToOne
    private BioentryEntity subjectBioEntry;
    
    @Column(name = "DISTANCE")
    private Integer distance;
    
    @JoinColumn(name = "TERM_ID", referencedColumnName = "TERM_ID", insertable = false, updatable = false)
    @ManyToOne
    private TermEntity term;

    public BioentryPathEntity() {
    }

    public BioentryPathEntity(BioentryPathEntityUK bioentryPathUK){
        this.id = bioentryPathUK;
    }
    
    public BioentryPathEntity(int objectBioentryId, int subjectBioentryId, int termId ){
        this.id = new BioentryPathEntityUK(objectBioentryId, subjectBioentryId, termId);
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (id != null ? id.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryPathEntity)) {
            return false;
        }
        BioentryPathEntity other = (BioentryPathEntity) object;
        if ((this.id == null && other.id != null) || (this.id != null && !this.id.equals(other.id))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryPath[id=" + id + "]";
    }

    public BioentryEntity getObjectBioEntry() {
        return objectBioEntry;
    }

    public void setObjectBioEntry(BioentryEntity objectBioEntry) {
        this.objectBioEntry = objectBioEntry;
    }

    public BioentryEntity getSubjectBioEntry() {
        return subjectBioEntry;
    }

    public void setSubjectBioEntry(BioentryEntity subjectBioEntry) {
        this.subjectBioEntry = subjectBioEntry;
    }

    public Integer getDistance() {
        return distance;
    }

    public void setDistance(Integer distance) {
        this.distance = distance;
    }

    public TermEntity getTerm() {
        return term;
    }

    public void setTerm(TermEntity term) {
        this.term = term;
    }

}
