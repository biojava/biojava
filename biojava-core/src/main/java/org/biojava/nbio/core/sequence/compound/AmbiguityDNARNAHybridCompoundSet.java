package org.biojava.nbio.core.sequence.compound;

/**
 * Ambiguity set for hybrid DNA/RNA sequences. Needed for some instances of synthetic nucleotide sequences present in protein structures from the PDB.
 * 
 * @author Jose Duarte
 *
 */
public class AmbiguityDNARNAHybridCompoundSet extends DNACompoundSet {

	private static class InitaliseOnDemand {
	    public static final AmbiguityDNARNAHybridCompoundSet INSTANCE = new AmbiguityDNARNAHybridCompoundSet();
	  }

	  public static AmbiguityDNARNAHybridCompoundSet getDNARNAHybridCompoundSet() {
	    return InitaliseOnDemand.INSTANCE;
	  }

	  public AmbiguityDNARNAHybridCompoundSet() {
	    super();

	    // this is the only one needed to make it a hybrid DNA/RNA. The rest are the usual DNA/RNA ambiguity letters
	    addNucleotideCompound("U", "A");
	    
	    
	    addNucleotideCompound("M", "K",
	        "A", "C");
	    addNucleotideCompound("R", "Y",
	        "A", "G");
	    addNucleotideCompound("W", "W",
	        "A", "T");
	    addNucleotideCompound("S", "S",
	        "C", "G");
	    addNucleotideCompound("Y", "R",
	        "C", "T");
	    addNucleotideCompound("K", "M",
	        "G", "T");
	    addNucleotideCompound("V", "B",
	        "A", "C", "G");
	    addNucleotideCompound("H", "D",
	        "A", "C", "T");
	    addNucleotideCompound("D", "H",
	        "A", "G", "T");
	    addNucleotideCompound("B", "V",
	        "C", "G", "T");
	    addNucleotideCompound("N", "N", "A", "C", "G", "T");

	    addNucleotideCompound("I", "I", "N", "A", "C", "G", "T");

	    
	    calculateIndirectAmbiguities();
	    
	  }
}
