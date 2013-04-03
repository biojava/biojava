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
 */
package org.biojava3.aaproperties.xml;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Set of proteinogenic amino acids.  Molecular weights are recorded in daltons (Da) as residues of a chain; monomers
 * outside of a chain would likely have an additional mass of 18.01524 Da contributed by an associated water molecule.
 *
 * Currently we have different symbols to handle inserts so not as clean as it should be
 *
 * @author Richard Holland
 * @author Scooter Willis
 * @author Mark Chapman
 */
public class CaseFreeAminoAcidCompoundSet implements CompoundSet<AminoAcidCompound> {

    private final Map<String, AminoAcidCompound> aminoAcidCompoundCache = new HashMap<String, AminoAcidCompound>();
    private final Map<AminoAcidCompound, Set<AminoAcidCompound>> equivalentsCache =
            new HashMap<AminoAcidCompound, Set<AminoAcidCompound>>();

    public CaseFreeAminoAcidCompoundSet() {
        aminoAcidCompoundCache.put("A", new AminoAcidCompound(null, "A", "Ala", "Alanine", 71.0788f));
        aminoAcidCompoundCache.put("R", new AminoAcidCompound(null, "R", "Arg", "Arginine", 156.1875f));
        aminoAcidCompoundCache.put("N", new AminoAcidCompound(null, "N", "Asn", "Asparagine", 114.1039f));
        aminoAcidCompoundCache.put("D", new AminoAcidCompound(null, "D", "Asp", "Aspartic acid", 115.0886f));
        aminoAcidCompoundCache.put("C", new AminoAcidCompound(null, "C", "Cys", "Cysteine", 103.1388f));
        aminoAcidCompoundCache.put("E", new AminoAcidCompound(null, "E", "Glu", "Glutamic acid", 129.1155f));
        aminoAcidCompoundCache.put("Q", new AminoAcidCompound(null, "Q", "Gln", "Glutamine", 128.1307f));
        aminoAcidCompoundCache.put("G", new AminoAcidCompound(null, "G", "Gly", "Glycine", 57.0519f));
        aminoAcidCompoundCache.put("H", new AminoAcidCompound(null, "H", "His", "Histidine", 137.1411f));
        aminoAcidCompoundCache.put("I", new AminoAcidCompound(null, "I", "Ile", "Isoleucine", 113.1594f));
        aminoAcidCompoundCache.put("L", new AminoAcidCompound(null, "L", "Leu", "Leucine", 113.1594f));
        aminoAcidCompoundCache.put("K", new AminoAcidCompound(null, "K", "Lys", "Lysine", 128.1741f));
        aminoAcidCompoundCache.put("M", new AminoAcidCompound(null, "M", "Met", "Methionine", 131.1986f));
        aminoAcidCompoundCache.put("F", new AminoAcidCompound(null, "F", "Phe", "Phenylalanine", 147.1766f));
        aminoAcidCompoundCache.put("P", new AminoAcidCompound(null, "P", "Pro", "Proline", 97.1167f));
        aminoAcidCompoundCache.put("S", new AminoAcidCompound(null, "S", "Ser", "Serine", 87.0782f));
        aminoAcidCompoundCache.put("T", new AminoAcidCompound(null, "T", "Thr", "Threonine", 101.1051f));
        aminoAcidCompoundCache.put("W", new AminoAcidCompound(null, "W", "Trp", "Tryptophan", 186.2132f));
        aminoAcidCompoundCache.put("Y", new AminoAcidCompound(null, "Y", "Tyr", "Tyrosine", 163.1760f));
        aminoAcidCompoundCache.put("V", new AminoAcidCompound(null, "V", "Val", "Valine", 99.1326f));
        aminoAcidCompoundCache.put("B", new AminoAcidCompound(null, "B", "Asx", "Asparagine or Aspartic acid", null));
        aminoAcidCompoundCache.put("Z", new AminoAcidCompound(null, "Z", "Glx", "Glutamine or Glutamic acid", null));
        aminoAcidCompoundCache.put("J", new AminoAcidCompound(null, "J", "Xle", "Leucine or Isoleucine", null));
        aminoAcidCompoundCache.put("X", new AminoAcidCompound(null, "X", "Xaa", "Unspecified", null));
        aminoAcidCompoundCache.put("-", new AminoAcidCompound(null, "-", "---", "Unspecified", null));
        aminoAcidCompoundCache.put(".", new AminoAcidCompound(null, ".", "...", "Unspecified", null));
        aminoAcidCompoundCache.put("_", new AminoAcidCompound(null, "_", "___", "Unspecified", null));
        aminoAcidCompoundCache.put("*", new AminoAcidCompound(null, "*", "***", "Stop", null));
        
        //Selenocysteine - this is encoded by UGA with the presence
        //of a SECIS element (SElenoCysteine Insertion Sequence) in the mRNA
        //and is a post-translation modification
        aminoAcidCompoundCache.put("U", new AminoAcidCompound(null, "U", "Sec", "Selenocysteine", 150.0388f));

        //Pyrrolysine is encoded by UAG in mRNA (normally Amber stop codon) which is translated to
        //this amino acid under the presence of pylT which creates an anti-codon CUA & pylS
        //which then does the actual conversion to Pyl.
        aminoAcidCompoundCache.put("O", new AminoAcidCompound(null, "O", "Pyl", "Pyrrolysine", 255.3172f));
        
        Map<String, AminoAcidCompound> lowerCaseSet = new HashMap<String, AminoAcidCompound>();
        for(String s:this.aminoAcidCompoundCache.keySet()){
        	lowerCaseSet.put(s.toLowerCase(), this.aminoAcidCompoundCache.get(s));
        }
        this.aminoAcidCompoundCache.putAll(lowerCaseSet);
    }

    public String getStringForCompound(AminoAcidCompound compound) {
        return compound.toString();
    }

    public AminoAcidCompound getCompoundForString(String string) {
        if (string.length() == 0) {
            return null;
        }
        if (string.length() > this.getMaxSingleCompoundStringLength()) {
            throw new IllegalArgumentException("String supplied ("+string+") is too long. Max is "+getMaxSingleCompoundStringLength());
        }
        return this.aminoAcidCompoundCache.get(string);
    }

    public int getMaxSingleCompoundStringLength() {
        return 1;
    }


    public boolean isCompoundStringLengthEqual() {
        return true;
    }

    private final static CaseFreeAminoAcidCompoundSet aminoAcidCompoundSet = new CaseFreeAminoAcidCompoundSet();

    public static CaseFreeAminoAcidCompoundSet getAminoAcidCompoundSet() {
        return aminoAcidCompoundSet;
    }

    public boolean compoundsEquivalent(AminoAcidCompound compoundOne, AminoAcidCompound compoundTwo) {
        Set<AminoAcidCompound> equivalents = getEquivalentCompounds(compoundOne);
        return (equivalents != null) && equivalents.contains(compoundTwo);
    }

    public Set<AminoAcidCompound> getEquivalentCompounds(AminoAcidCompound compound) {
        if (equivalentsCache.isEmpty()) {
            // most compounds are equivalent to themselves alone
            for (AminoAcidCompound c : aminoAcidCompoundCache.values()) {
                equivalentsCache.put(c, Collections.singleton(c));
            }
            // ambiguous Asparagine or Aspartic acid
            addAmbiguousEquivalents("N", "D", "B");
            // ambiguous Glutamine or Glutamic acid
            addAmbiguousEquivalents("E", "Q", "Z");
            // ambiguous Leucine or Isoleucine
            addAmbiguousEquivalents("I", "L", "J");
            // ambiguous gaps
            AminoAcidCompound gap1, gap2, gap3;
            Set<AminoAcidCompound> gaps = new HashSet<AminoAcidCompound>();
            gaps.add(gap1 = aminoAcidCompoundCache.get("-"));
            gaps.add(gap2 = aminoAcidCompoundCache.get("."));
            gaps.add(gap3 = aminoAcidCompoundCache.get("_"));
            equivalentsCache.put(gap1, gaps);
            equivalentsCache.put(gap2, gaps);
            equivalentsCache.put(gap3, gaps);
            // X is never equivalent, even to itself
            equivalentsCache.put(aminoAcidCompoundCache.get("X"), new HashSet<AminoAcidCompound>());
        }
        return equivalentsCache.get(compound);
    }

    // helper method to initialize the equivalent sets for 2 amino acid compounds and their ambiguity compound
    private void addAmbiguousEquivalents(String one, String two, String either) {
        Set<AminoAcidCompound> equivalents;
        AminoAcidCompound cOne, cTwo, cEither;

        equivalents = new HashSet<AminoAcidCompound>();
        equivalents.add(cOne = aminoAcidCompoundCache.get(one));
        equivalents.add(cTwo = aminoAcidCompoundCache.get(two));
        equivalents.add(cEither = aminoAcidCompoundCache.get(either));
        equivalentsCache.put(cEither, equivalents);

        equivalents = new HashSet<AminoAcidCompound>();
        equivalents.add(cOne);
        equivalents.add(cEither);
        equivalentsCache.put(cOne, equivalents);

        equivalents = new HashSet<AminoAcidCompound>();
        equivalents.add(cTwo);
        equivalents.add(cEither);
        equivalentsCache.put(cTwo, equivalents);
    }

    public boolean hasCompound(AminoAcidCompound compound) {
        return aminoAcidCompoundCache.containsValue(compound);
    }

    // TODO throwing an error seems unnecessary, should this return a boolean instead? maybe rename to isValidSequence?
    public void verifySequence(Sequence<AminoAcidCompound> sequence) throws CompoundNotFoundError {
        for (AminoAcidCompound compound : sequence) {
            if (!hasCompound(compound)) {
                throw new CompoundNotFoundError("Compound (" + compound + ") not found in AminoAcidCompoundSet.");
            }
        }
    }

    public List<AminoAcidCompound> getAllCompounds() {
        return new ArrayList<AminoAcidCompound>(aminoAcidCompoundCache.values());
    }


    public boolean isComplementable() {
        return false;
    }
}
