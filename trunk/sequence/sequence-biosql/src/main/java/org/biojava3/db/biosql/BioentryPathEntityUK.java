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

import javax.persistence.Column;
import javax.persistence.Embeddable;

/**
 * The embedded unique key for BioentryPathEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class BioentryPathEntityUK implements java.io.Serializable {
    private static final long serialVersionUID = -1434846983188237697L;
    @Column(name = "TERM_ID", nullable = false)
    private int termId;
    @Column(name = "OBJECT_BIOENTRY_ID", nullable = false)
    private int objectBioentryId;
    @Column(name = "SUBJECT_BIOENTRY_ID", nullable = false)
    private int subjectBioentryId;
    
    public BioentryPathEntityUK(){}
    
    public BioentryPathEntityUK(int objectBioentryId, int subjectBioentryId, int termId){
        this.objectBioentryId = objectBioentryId;
        this.subjectBioentryId = subjectBioentryId;
        this.termId = termId;
    }

    public int getTermId() {
        return termId;
    }

    public void setTermId(int termId) {
        this.termId = termId;
    }

    public int getObjectBioentryId() {
        return objectBioentryId;
    }

    public void setObjectBioentryId(int objectBioentryId) {
        this.objectBioentryId = objectBioentryId;
    }

    public int getSubjectBioentryId() {
        return subjectBioentryId;
    }

    public void setSubjectBioentryId(int subjectBioentryId) {
        this.subjectBioentryId = subjectBioentryId;
    }
    
    @Override
    public int hashCode() {
        int hash = 0;
        hash += (int) objectBioentryId;
        hash += (int) termId;
        hash += (int) subjectBioentryId;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BioentryPathEntityUK)) {
            return false;
        }
        BioentryPathEntityUK other = (BioentryPathEntityUK) object;
        if (this.objectBioentryId != other.objectBioentryId) {
            return false;
        }
        if (this.termId != other.termId) {
            return false;
        }
        if (this.subjectBioentryId != other.subjectBioentryId) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryPathUK[objectBioentryId=" + 
                objectBioentryId + ", termId=" + termId + ", subjectBioentryId=" + subjectBioentryId + "]";
    }

}
