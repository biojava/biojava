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
 * Created on 01-21-2010
 *
 * @author Richard Holland
 * @auther Scooter Willis
 *
 */
package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.NucleotideCompoundInterface;

public class NucleotideCompound<C extends Compound> implements NucleotideCompoundInterface {

    private String base;
    private CompoundSet<NucleotideCompound> compoundSet;
    private String complementStr;

    public NucleotideCompound(String base, CompoundSet<NucleotideCompound> compoundSet, String complementStr) {
        this.base = base;
        this.compoundSet = compoundSet;
        this.complementStr = complementStr;
    }

    public NucleotideCompound getComplement() {
        return this.compoundSet.getCompoundForString(this.complementStr);
    }

    public String toString() {
        return base.toString();
    }

    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof NucleotideCompound)) {
            return false;
        }
        NucleotideCompound them = (NucleotideCompound) obj;
        return this.base.equals(them.base);
    }

    public int hashCode() {
        return this.base.hashCode();
    }

    public boolean equalsIgnoreCase(Compound compound) {
        if (compound == null) {
            return false;
        }
        if (!(compound instanceof NucleotideCompound)) {
            return false;
        }
        NucleotideCompound them = (NucleotideCompound) compound;
        return this.base.toString().equalsIgnoreCase(them.base.toString());
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

    public void setComplement(NucleotideCompoundInterface nucleotideCompound) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
