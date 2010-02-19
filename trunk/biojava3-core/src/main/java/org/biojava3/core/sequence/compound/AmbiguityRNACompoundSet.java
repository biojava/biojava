package org.biojava3.core.sequence.compound;

/**
 *
 * @author Andy Yates
 */
public class AmbiguityRNACompoundSet extends RNACompoundSet {

  private static class InitaliseOnDemand {
    public static final AmbiguityRNACompoundSet INSTANCE = new AmbiguityRNACompoundSet();
  }

  public static AmbiguityRNACompoundSet getDNACompoundSet() {
    return InitaliseOnDemand.INSTANCE;
  }

  public AmbiguityRNACompoundSet() {

    addNucleotideCompound("A", "U");
    addNucleotideCompound("U", "A");
    addNucleotideCompound("G", "C");
    addNucleotideCompound("C", "G");

    addNucleotideCompound("M", "K",
        "A", "C");
    addNucleotideCompound("R", "Y",
        "A", "G");
    addNucleotideCompound("W", "W",
        "A", "U");
    addNucleotideCompound("S", "S",
        "C", "G");
    addNucleotideCompound("Y", "R",
        "C", "U");
    addNucleotideCompound("K", "M",
        "G", "U");
    addNucleotideCompound("V", "B",
        "A", "C", "G");
    addNucleotideCompound("H", "D",
        "A", "C", "U");
    addNucleotideCompound("D", "H",
        "A", "G", "U");
    addNucleotideCompound("B", "V",
        "C", "G", "U");
    addNucleotideCompound("N", "N",
        "A", "C", "G", "U", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B");

    calculateIndirectAmbiguities();
  }

}
