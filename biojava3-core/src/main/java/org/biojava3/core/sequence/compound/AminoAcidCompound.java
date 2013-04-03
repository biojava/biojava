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

import org.biojava3.core.sequence.template.AbstractCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;

/**
 * Used to describe an Amino Acid.
 * @author Richard Holland
 * @author Scooter Willis
 * @author Andy Yates
 */
public class AminoAcidCompound extends AbstractCompound {

  private final AminoAcidCompoundSet compoundSet;

  public AminoAcidCompound(AminoAcidCompoundSet compoundSet, String shortName,
      String longName, String description, Float molecularWeight) {
    super(shortName);
    setShortName(shortName);
    setLongName(longName);
    setDescription(description);
    setMolecularWeight(molecularWeight);
    this.compoundSet = compoundSet;
  }

  // TODO need to allow for modified name; that's not equality though is it?
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (!(obj instanceof AminoAcidCompound)) {
      return false;
    }
    AminoAcidCompound them = (AminoAcidCompound) obj;
    if (toString().equals(them.toString())) {
      return true;
    }
    return getLongName().equals(them.getLongName());

  }

  public int hashCode() {
    return toString().hashCode();
  }

  public boolean equalsIgnoreCase(Compound compound) {
    if (compound == null) {
      return false;
    }
    if (!(compound instanceof AminoAcidCompound)) {
      return false;
    }
    AminoAcidCompound them = (AminoAcidCompound) compound;
    if (toString().equalsIgnoreCase(them.toString())) {
      return true;
    }
    return getLongName().equalsIgnoreCase(them.getLongName());
  }

  public CompoundSet<AminoAcidCompound> getCompoundSet() {
    return compoundSet;
  }
}
