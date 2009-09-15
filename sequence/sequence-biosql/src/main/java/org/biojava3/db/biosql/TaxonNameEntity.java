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

/**
 * Entity for TaxonName
 * @author Mark Schreiber
 * @since v3
 */
@Entity
public class TaxonNameEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @EmbeddedId
    private TaxonNameEntityUK taxonNameUK;
    
    @JoinColumn(name = "TAXON_ID", referencedColumnName = "TAXON_ID", insertable = false, updatable = false)
    @ManyToOne
    private TaxonEntity taxon;
    @Column(name = "NAME")
    private String name;
    @Column(name = "NAME_CLASS")
    private String nameClass;
    
    public TaxonNameEntity(){}
    public TaxonNameEntity(int taxonId, String name, String nameClass){
        this.taxonNameUK = new TaxonNameEntityUK(taxonId, name, nameClass);
    }

    public void setTaxonNameUK(TaxonNameEntityUK id) {
        this.taxonNameUK = id;
    }

    public TaxonNameEntityUK getTaxonNameUK() {
        return taxonNameUK;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (taxonNameUK != null ? taxonNameUK.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the taxonNameUK fields are not set
        if (!(object instanceof TaxonNameEntity)) {
            return false;
        }
        TaxonNameEntity other = (TaxonNameEntity) object;
        if ((this.taxonNameUK == null && other.taxonNameUK != null) || (this.taxonNameUK != null && !this.taxonNameUK.equals(other.taxonNameUK))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.TaxonNameEntity [id=" + taxonNameUK + "]";
    }

    public TaxonEntity getTaxon() {
        return taxon;
    }

    public void setTaxon(TaxonEntity taxon) {
        this.taxon = taxon;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getNameClass() {
        return nameClass;
    }

    public void setNameClass(String nameClass) {
        this.nameClass = nameClass;
    }

}
