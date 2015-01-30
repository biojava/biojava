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
package org.biojava.nbio.core.sequence.compound;


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
