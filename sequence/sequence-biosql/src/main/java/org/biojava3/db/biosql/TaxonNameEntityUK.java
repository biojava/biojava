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
 * Unique key for TaxonNameEntity
 * @author Mark Schreiber
 * @since v3
 */
@Embeddable
public class TaxonNameEntityUK implements java.io.Serializable {
    private static final long serialVersionUID = 7410475607776666220L;
    @Column(name="taxon_id", nullable=false)
    private int taxonId;
    @Column(name="name", nullable=false, updatable=false)
    private String name;
    @Column(name="nameClass", nullable=false, updatable=false)
    private String nameClass;
    
    public TaxonNameEntityUK(){
        
    }
    public TaxonNameEntityUK(int taxonId, String name, String nameClass){
        this.taxonId = taxonId;
        this.name = name;
        this.nameClass = nameClass;
    }

    public int getTaxonId() {
        return taxonId;
    }

    public void setTaxonId(int taxonId) {
        this.taxonId = taxonId;
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
    
    @Override
    public String toString(){
        return "org.biosql.entity.TaxonNameEntityUK[taxonId=" + taxonId + ", name=" + name +", nameClass= "+ nameClass+"]";
    }
    
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final TaxonNameEntityUK other = (TaxonNameEntityUK) obj;
        if (this.taxonId != other.taxonId) {
            return false;
        }
        if (this.name.equals(other.name) && (this.name == null || !this.name.equals(other.name))) {
            return false;
        }
        if (this.nameClass.equals(other.nameClass) && (this.nameClass == null || !this.nameClass.equals(other.nameClass))) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 79 * hash + this.taxonId;
        hash = 79 * hash + (this.name != null ? this.name.hashCode() : 0);
        hash = 79 * hash + (this.nameClass != null ? this.nameClass.hashCode() : 0);
        return hash;
    }
}
