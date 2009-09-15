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
import javax.persistence.Lob;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.OneToOne;
import javax.persistence.Table;

/**
 * The persistent entity for a BiosequenceEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "BIOSEQUENCE")
@NamedQueries({
    @NamedQuery(name = "BiosequenceEntity.findByBioentryId", query = "SELECT b FROM BiosequenceEntity b WHERE b.bioentryId = :bioentryId")})
    public class BiosequenceEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "BIOENTRY_ID", nullable = false)
    private Integer bioentryId;
    @Column(name = "VERSION")
    private Integer version;
    @Column(name = "LENGTH") //length is a JPA reserved word. This might cause issues
    private Integer length;
    @Column(name = "ALPHABET")
    private String alphabet;
    @Lob
    @Column(name = "SEQ")
    private String seq;
    @JoinColumn(name = "BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID", insertable = false, updatable = false)
    @OneToOne
    private BioentryEntity bioentry;

    public BiosequenceEntity() {
    }

    public BiosequenceEntity(Integer bioentryId) {
        this.bioentryId = bioentryId;
    }

    public Integer getBioentryId() {
        return bioentryId;
    }

    public void setBioentryId(Integer bioentryId) {
        this.bioentryId = bioentryId;
    }

    public Integer getVersion() {
        return version;
    }

    public void setVersion(Integer version) {
        this.version = version;
    }

    public Integer getLength() {
        return length;
    }

    public void setLength(Integer length) {
        this.length = length;
    }

    public String getAlphabet() {
        return alphabet;
    }

    public void setAlphabet(String alphabet) {
        this.alphabet = alphabet;
    }

    public String getSeq() {
        return seq;
    }

    public void setSeq(String seq) {
        this.seq = seq;
    }

    public BioentryEntity getBioentry() {
        return bioentry;
    }

    public void setBioentry(BioentryEntity bioentry) {
        this.bioentry = bioentry;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (bioentryId != null ? bioentryId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof BiosequenceEntity)) {
            return false;
        }
        BiosequenceEntity other = (BiosequenceEntity) object;
        if ((this.bioentryId == null && other.bioentryId != null) || (this.bioentryId != null && !this.bioentryId.equals(other.bioentryId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BiosequenceEntity[bioentryId=" + bioentryId + "]";
    }

}
