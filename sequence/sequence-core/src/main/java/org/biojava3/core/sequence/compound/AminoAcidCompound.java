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
 * Created on DATE
 *
 */

package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.Compound;

/* @author Richard Holland
 * @author Scooter Willis
 */

public class AminoAcidCompound implements Compound {

    
    private AminoAcidCompoundSet compoundSet;
    

    public AminoAcidCompound(AminoAcidCompoundSet compoundSet,String shortName,String longName,String description,Float molecularWeight ) {
        setShortName(shortName);
        setLongName(longName);
        setDescription(description);
        setMolecularWeight(molecularWeight);
        this.compoundSet = compoundSet;
        
    }



    public String toString() {
        return shortName.toString();
    }
//TODO need to allow for modified name
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof AminoAcidCompound)) {
            return false;
        }
        AminoAcidCompound them = (AminoAcidCompound) obj;
        if(shortName.equals(them.shortName)){
            return true;
        }
        return longName.equals(them.longName);

    }

//    public int hashCode() {
//        return this.base.hashCode();
//    }

    public boolean equalsIgnoreCase(Compound compound) {
        if (compound == null) {
            return false;
        }
        if (!(compound instanceof AminoAcidCompound)) {
            return false;
        }
        AminoAcidCompound them = (AminoAcidCompound) compound;
        if(shortName.toString().equalsIgnoreCase(them.shortName.toString())){
            return true;
        }
        return longName.toString().equalsIgnoreCase(them.shortName.toString());
    }

    private String shortName = null;
    private String longName = null;
    private String description = null;
    private Float molecularWeight = null;

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getShortName() {
        return shortName;
    }

    public void setShortName(String shortName) {
        this.shortName = shortName;
    }

    public String getLongName() {
        return longName;
    }

    public void setLongName(String longName) {
        this.longName = longName;
    }

    public Float getMolecularWeight() {
        return molecularWeight;
    }

    public void setMolecularWeight(Float molecularWeight) {
        this.molecularWeight = molecularWeight;
    }
}
