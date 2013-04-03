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

import static java.util.Arrays.asList;
import static java.util.Collections.unmodifiableSet;

import java.util.HashSet;
import java.util.Set;

import org.biojava3.core.sequence.template.AbstractCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ComplementCompound;

/**
 *
 * @author Scooter Willis
 * @author Andy Yates
 */
public class NucleotideCompound extends AbstractCompound implements ComplementCompound {

    private final CompoundSet<NucleotideCompound> compoundSet;
    private final String complementStr;
    private final Set<NucleotideCompound> constituents;

    public NucleotideCompound(String base, CompoundSet<NucleotideCompound> compoundSet, String complementStr) {
      super(base);
      this.compoundSet = compoundSet;
      this.complementStr = complementStr;
      this.constituents = unmodifiableSet(new HashSet<NucleotideCompound>(asList(this)));
    }

    public NucleotideCompound(String base, CompoundSet<NucleotideCompound> compoundSet, String complementStr, NucleotideCompound[] constituents) {
        super(base);
        this.compoundSet = compoundSet;
        this.complementStr = complementStr;
        this.constituents = unmodifiableSet(new HashSet<NucleotideCompound>(asList(constituents)));
    }

    @Override
    public String getShortName() {
      return getBase();
    }

    public ComplementCompound getComplement() {
        return compoundSet.getCompoundForString(complementStr);
    }

    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof NucleotideCompound)) {
            return false;
        }
        NucleotideCompound them = (NucleotideCompound) obj;
        return toString().equals(them.toString());
    }

    public int hashCode() {
        return toString().hashCode();
    }

    public boolean equalsIgnoreCase(Compound compound) {
        if (compound == null) {
            return false;
        }
        if (!(compound instanceof NucleotideCompound)) {
            return false;
        }
        NucleotideCompound them = (NucleotideCompound) compound;
        return toString().equalsIgnoreCase(them.toString());
    }

    public Set<NucleotideCompound> getConstituents() {
      return constituents;
    }
    
    /**@deprecated replaced with {@link #getConstituents()} due to typographical error */
    public Set<NucleotideCompound> getConsituents() {
    	return getConstituents();
    }

    public boolean isAmbiguous() {
      return !constituents.isEmpty();
    }
}
