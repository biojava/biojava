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
 *
 */
package org.biojava.nbio.structure.basepairs;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.contact.Pair;

import javax.vecmath.Matrix4d;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * This class allows for finding inter-strand base pairs that are not necessarily canonical Watson-Crick pairs.
 * The implementation of findPair is different than that of the base class.  This class does not consider intra-strand
 * base pairing and for that, the TertiaryBasePairParameters class should be used.
 * @author Luke Czapla
 * @since 5.0.0
 *
 */
public class MismatchedBasePairParameters extends BasePairParameters implements Serializable {

    // These are the criteria used to select proper base pairs.
    protected static double MaxStagger = 2.0, MaxShear = 5.0, MaxStretch = 5.0,
            MaxPropeller = 60.0;

    /**
     * This constructor is used to create the TertiaryBasePairParameters object.  The parent constructors are valid
     * as well, but for this class, it makes the most sense to specify the exact parameters for the analysis.
     * @param structure The Structure to analyze
     * @param RNA Whether to analyze RNA (if false, it will analyze DNA)
     * @param removeDups Whether to remove duplicate sequences (useful for RCSB data with redundant units).
     * @param canonical Whether to only consider canonical Watson-Crick base pairs.  If false, any pairing will be identified
     *                  as long it falls below the maximum values of stagger, shear, and stretch.
     */
    public MismatchedBasePairParameters(Structure structure, boolean RNA, boolean removeDups, boolean canonical) {

        super(structure, RNA, removeDups, canonical);

    }

    /**
     * This is an implementation for finding non-canonical base pairs when there may be missing or overhanging bases.
     * @param chains The list of chains already found to be nucleic acids.
     * @return The list of the atom groups (residues) that are pairs, as a Pair of nucleic acid Groups.
     */
    @Override
    public List<Pair<Group>> findPairs(List<Chain> chains) {
        List<Pair<Group>> result = new ArrayList<>();
        boolean lastFoundPair = false;
        for (int i = 0; i < chains.size(); i++) {
            Chain c = chains.get(i);
            String sequence = c.getAtomSequence();
            for (int m = 0; m < sequence.length(); m++) {
                boolean foundPair = false;
                Integer type1, type2;
                for (int j = i + 1; j < chains.size() && !foundPair; j++) {
                    Chain c2 = chains.get(j);
                    if (j > i+1 && c.getAtomSequence().equals(c2.getAtomSequence()) && nonredundant) continue;
                    String sequence2 = c2.getAtomSequence();
                    for (int k = c2.getAtomSequence().length() - 1; k >= 0 && !foundPair; k--) {
                        if (canonical && !BasePairParameters.match(sequence.charAt(m), sequence2.charAt(k), useRNA)) continue;
                        Group g1 = c.getAtomGroup(m);
                        Group g2 = c2.getAtomGroup(k);
                        type1 = BASE_MAP.get(g1.getPDBName());
                        type2 = BASE_MAP.get(g2.getPDBName());
                        if (type1 == null || type2 == null) continue;
                        Atom a1 = g1.getAtom("C1'");
                        Atom a2 = g2.getAtom("C1'");
                        if (a1 == null || a2 == null) continue;
                        // C1'-C1' distance is one useful criteria
                        if (Math.abs(a1.getCoordsAsPoint3d().distance(a2.getCoordsAsPoint3d()) - 10.0) > 4.0) continue;
                        Pair<Group> ga = new Pair<>(g1, g2);
                        Matrix4d data = basePairReferenceFrame(ga);
                        // if the stagger is greater than 2 Å, it's not really paired.
                        if (Math.abs(pairParameters[5]) > MaxStagger) continue;
                        // similarly, extreme shear and stretch is not a good base pair
                        if (Math.abs(pairParameters[3]) > MaxShear) continue;
                        if (Math.abs(pairParameters[4]) > MaxStretch) continue;

                        // if the propeller is ridiculous it's also not that good of a pair.
                        if (Math.abs(pairParameters[1]) > MaxPropeller) {
                            continue;
                        }
                        result.add(ga);
                        pairingNames.add(useRNA ? BASE_LIST_RNA[type1] + BASE_LIST_RNA[type2] : BASE_LIST_DNA[type1] + BASE_LIST_DNA[type2]);
                        foundPair = true;
                    }
                    if (!foundPair && lastFoundPair) {
                        if (pairSequence.length() > 0 && pairSequence.charAt(pairSequence.length() - 1) != ' ')
                            pairSequence += ' ';
                    }
                    if (foundPair) pairSequence += (c.getAtomSequence().charAt(i));
                    lastFoundPair = foundPair;
                }
            }
        }
        return result;
    }

    public static double getMaxStagger() {
        return MaxStagger;
    }

    public static void setMaxStagger(double maxStagger) {
        MaxStagger = maxStagger;
    }

    public static double getMaxShear() {
        return MaxShear;
    }

    public static void setMaxShear(double maxShear) {
        MaxShear = maxShear;
    }

    public static double getMaxStretch() {
        return MaxStretch;
    }

    public static void setMaxStretch(double maxStretch) {
        MaxStretch = maxStretch;
    }

    public static double getMaxPropeller() {
        return MaxPropeller;
    }

    public static void setMaxPropeller(double maxPropeller) {
        MaxPropeller = maxPropeller;
    }
}
