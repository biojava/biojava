package org.biojava3.core.sequence.compound;

import org.biojava3.core.sequence.template.AbstractNucleotideCompoundSet;

/**
 *
 * @author Andy Yates
 */
public class RNACompoundSet extends AbstractNucleotideCompoundSet<NucleotideCompound> {

  private static class InitaliseOnDemand {
    public static final RNACompoundSet INSTANCE = new RNACompoundSet();
  }
  public static RNACompoundSet getRNACompoundSet() {
    return InitaliseOnDemand.INSTANCE;
  }

  public RNACompoundSet() {
    addNucleotideCompound("A", "U");
    addNucleotideCompound("U", "A");
    addNucleotideCompound("G", "C");
    addNucleotideCompound("C", "G");
    addNucleotideCompound("N", "N");
    addNucleotideCompound("-", "-");
  }

  public NucleotideCompound newNucleotideCompound(String base, String complement, String... equivalents) {
    return new NucleotideCompound(base, this, complement);
  }
}