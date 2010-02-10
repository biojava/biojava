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
 * The entity bean that backs a BioentryEntity implementation
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "BIOENTRY")
@NamedQueries({
    @NamedQuery(name = "BioentryEntity.findByBioentryId", query = "SELECT b FROM BioentryEntity b WHERE b.bioentryId = :bioentryId"), 
    @NamedQuery(name = "BioentryEntity.findByName", query = "SELECT b FROM BioentryEntity b WHERE b.name = :name"), 
    @NamedQuery(name = "BioentryEntity.findByAccession", query = "SELECT b FROM BioentryEntity b WHERE b.accession = :accession"), 
    @NamedQuery(name = "BioentryEntity.findByIdentifier", query = "SELECT b FROM BioentryEntity b WHERE b.identifier = :identifier"), 
    @NamedQuery(name = "BioentryEntity.findByDivision", query = "SELECT b FROM BioentryEntity b WHERE b.division = :division"), 
    @NamedQuery(name = "BioentryEntity.findByDescription", query = "SELECT b FROM BioentryEntity b WHERE b.description = :description"), 
    @NamedQuery(name = "BioentryEntity.findByVersion", query = "SELECT b FROM BioentryEntity b WHERE b.version = :version")})
public class BioentryEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "BIOENTRY_ID", nullable = false)
    private Integer bioentryId;
    @Column(name = "NAME", nullable = false)
    private String name;
    @Column(name = "ACCESSION", nullable = false)
    private String accession;
    @Column(name = "IDENTIFIER", nullable = false)
    private String identifier;
    @Column(name = "DIVISION")
    private String division;
    @Column(name = "DESCRIPTION")
    private String description;
    @Column(name = "VERSION", nullable = false)
    private int version;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "bioentry")
    private Collection<BioentryDbxrefEntity> bioentryDbxrefCollection;
    @JoinColumn(name = "BIODATABASE_ID", referencedColumnName = "BIODATABASE_ID")
    @ManyToOne
    private BiodatabaseEntity biodatabaseId;
    @JoinColumn(name = "TAXON_ID", referencedColumnName = "TAXON_ID")
    @ManyToOne
    private TaxonEntity taxonId;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "bioentryId")
    private Collection<CommentEntity> commentCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "bioentryId")
    private Collection<SeqfeatureEntity> seqfeatureCollection;
    @OneToOne(cascade = CascadeType.ALL, mappedBy = "bioentry")
    private BiosequenceEntity biosequence;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "subjectBioentryId")
    private Collection<BioentryRelationshipEntity> bioentryRelationshipCollection;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "objectBioentryId")
    private Collection<BioentryRelationshipEntity> bioentryRelationshipCollection1;
    @OneToMany(cascade = CascadeType.ALL, mappedBy = "bioentry")
    private Collection<BioentryReferenceEntity> bioentryReferenceCollection;

    public BioentryEntity() {
    }

    public BioentryEntity(Integer bioentryId) {
        this.bioentryId = bioentryId;
    }

    public BioentryEntity(Integer bioentryId, String name, String accession, String identifier, int version) {
        this.bioentryId = bioentryId;
        this.name = name;
        this.accession = accession;
        this.identifier = identifier;
        this.version = version;
    }

    public Integer getBioentryId() {
        return bioentryId;
    }

    public void setBioentryId(Integer bioentryId) {
        this.bioentryId = bioentryId;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getIdentifier() {
        return identifier;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    public String getDivision() {
        return division;
    }

    public void setDivision(String division) {
        this.division = division;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public int getVersion() {
        return version;
    }

    public void setVersion(int version) {
        this.version = version;
    }

    public Collection<BioentryDbxrefEntity> getBioentryDbxrefCollection() {
        return bioentryDbxrefCollection;
    }

    public void setBioentryDbxrefCollection(Collection<BioentryDbxrefEntity> bioentryDbxrefCollection) {
        this.bioentryDbxrefCollection = bioentryDbxrefCollection;
    }

    public BiodatabaseEntity getBiodatabaseId() {
        return biodatabaseId;
    }

    public void setBiodatabaseId(BiodatabaseEntity biodatabaseId) {
        this.biodatabaseId = biodatabaseId;
    }

    public TaxonEntity getTaxonId() {
        return taxonId;
    }

    public void setTaxonId(TaxonEntity taxonId) {
        this.taxonId = taxonId;
    }

    public Collection<CommentEntity> getCommentCollection() {
        return commentCollection;
    }

    public void setCommentCollection(Collection<CommentEntity> commentCollection) {
        this.commentCollection = commentCollection;
    }

    public Collection<SeqfeatureEntity> getSeqfeatureCollection() {
        return seqfeatureCollection;
    }

    public void setSeqfeatureCollection(Collection<SeqfeatureEntity> seqfeatureCollection) {
        this.seqfeatureCollection = seqfeatureCollection;
    }

    public BiosequenceEntity getBiosequence() {
        return biosequence;
    }

    public void setBiosequence(BiosequenceEntity biosequence) {
        this.biosequence = biosequence;
    }

    public Collection<BioentryRelationshipEntity> getBioentryRelationshipCollection() {
        return bioentryRelationshipCollection;
    }

    public void setBioentryRelationshipCollection(Collection<BioentryRelationshipEntity> bioentryRelationshipCollection) {
        this.bioentryRelationshipCollection = bioentryRelationshipCollection;
    }

    public Collection<BioentryRelationshipEntity> getBioentryRelationshipCollection1() {
        return bioentryRelationshipCollection1;
    }

    public void setBioentryRelationshipCollection1(Collection<BioentryRelationshipEntity> bioentryRelationshipCollection1) {
        this.bioentryRelationshipCollection1 = bioentryRelationshipCollection1;
    }

    public Collection<BioentryReferenceEntity> getBioentryReferenceCollection() {
        return bioentryReferenceCollection;
    }

    public void setBioentryReferenceCollection(Collection<BioentryReferenceEntity> bioentryReferenceCollection) {
        this.bioentryReferenceCollection = bioentryReferenceCollection;
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
        if (!(object instanceof BioentryEntity)) {
            return false;
        }
        BioentryEntity other = (BioentryEntity) object;
        if ((this.bioentryId == null && other.bioentryId != null) || (this.bioentryId != null && !this.bioentryId.equals(other.bioentryId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.BioentryEntity[bioentryId=" + bioentryId + "]";
    }

}
