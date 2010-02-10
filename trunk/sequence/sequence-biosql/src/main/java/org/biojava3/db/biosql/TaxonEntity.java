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
 * Entity for TaxonEntity
 * @author Mark Schreiber
 * @since v3
 */
@Entity
@Table(name = "TAXON")
@NamedQueries({
    @NamedQuery(name = "TaxonEntity.findByTaxonId", query = "SELECT t FROM TaxonEntity t WHERE t.taxonId = :taxonId"), 
    @NamedQuery(name = "TaxonEntity.findByNcbiTaxonId", query = "SELECT t FROM TaxonEntity t WHERE t.ncbiTaxonId = :ncbiTaxonId"), 
    @NamedQuery(name = "TaxonEntity.findByParentTaxonId", query = "SELECT t FROM TaxonEntity t WHERE t.parentTaxonId = :parentTaxonId"), 
    @NamedQuery(name = "TaxonEntity.findByNodeRank", query = "SELECT t FROM TaxonEntity t WHERE t.nodeRank = :nodeRank"), 
    @NamedQuery(name = "TaxonEntity.findByGeneticCode", query = "SELECT t FROM TaxonEntity t WHERE t.geneticCode = :geneticCode"), 
    @NamedQuery(name = "TaxonEntity.findByMitoGeneticCode", query = "SELECT t FROM TaxonEntity t WHERE t.mitoGeneticCode = :mitoGeneticCode"), 
    @NamedQuery(name = "TaxonEntity.findByLeftValue", query = "SELECT t FROM TaxonEntity t WHERE t.leftValue = :leftValue"), 
    @NamedQuery(name = "TaxonEntity.findByRightValue", query = "SELECT t FROM TaxonEntity t WHERE t.rightValue = :rightValue")})
public class TaxonEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "TAXON_ID", nullable = false)
    private Integer taxonId;
    @Column(name = "NCBI_TAXON_ID", nullable = false)
    private int ncbiTaxonId;
    @Column(name = "PARENT_TAXON_ID")
    private Integer parentTaxonId;
    @Column(name = "NODE_RANK")
    private String nodeRank;
    @Column(name = "GENETIC_CODE")
    private Short geneticCode;
    @Column(name = "MITO_GENETIC_CODE")
    private Short mitoGeneticCode;
    @Column(name = "LEFT_VALUE", nullable = false)
    private int leftValue;
    @Column(name = "RIGHT_VALUE", nullable = false)
    private int rightValue;
    @OneToMany(mappedBy = "taxonId")
    private Collection<BioentryEntity> bioentryCollection;

    public TaxonEntity() {
    }

    public TaxonEntity(Integer taxonId) {
        this.taxonId = taxonId;
    }

    public TaxonEntity(Integer taxonId, int ncbiTaxonId, int leftValue, int rightValue) {
        this.taxonId = taxonId;
        this.ncbiTaxonId = ncbiTaxonId;
        this.leftValue = leftValue;
        this.rightValue = rightValue;
    }

    public Integer getTaxonId() {
        return taxonId;
    }

    public void setTaxonId(Integer taxonId) {
        this.taxonId = taxonId;
    }

    public int getNcbiTaxonId() {
        return ncbiTaxonId;
    }

    public void setNcbiTaxonId(int ncbiTaxonId) {
        this.ncbiTaxonId = ncbiTaxonId;
    }

    public Integer getParentTaxonId() {
        return parentTaxonId;
    }

    public void setParentTaxonId(Integer parentTaxonId) {
        this.parentTaxonId = parentTaxonId;
    }

    public String getNodeRank() {
        return nodeRank;
    }

    public void setNodeRank(String nodeRank) {
        this.nodeRank = nodeRank;
    }

    public Short getGeneticCode() {
        return geneticCode;
    }

    public void setGeneticCode(Short geneticCode) {
        this.geneticCode = geneticCode;
    }

    public Short getMitoGeneticCode() {
        return mitoGeneticCode;
    }

    public void setMitoGeneticCode(Short mitoGeneticCode) {
        this.mitoGeneticCode = mitoGeneticCode;
    }

    public int getLeftValue() {
        return leftValue;
    }

    public void setLeftValue(int leftValue) {
        this.leftValue = leftValue;
    }

    public int getRightValue() {
        return rightValue;
    }

    public void setRightValue(int rightValue) {
        this.rightValue = rightValue;
    }

    public Collection<BioentryEntity> getBioentryCollection() {
        return bioentryCollection;
    }

    public void setBioentryCollection(Collection<BioentryEntity> bioentryCollection) {
        this.bioentryCollection = bioentryCollection;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (taxonId != null ? taxonId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof TaxonEntity)) {
            return false;
        }
        TaxonEntity other = (TaxonEntity) object;
        if ((this.taxonId == null && other.taxonId != null) || (this.taxonId != null && !this.taxonId.equals(other.taxonId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TaxonEntity[taxonId=" + taxonId + "]";
    }

}
