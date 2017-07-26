package org.biojava.nbio.structure.basepairs;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.contact.Pair;

import javax.vecmath.Matrix4d;
import java.util.ArrayList;
import java.util.List;

/**
 * Contributed to BioJava under it's LGPL
 * This class also finds the base pairing and base-pair step parameters but has a broader definition
 * of a base pair so that non-canonical-WC base pairs will be detected and reported.  This is useful
 * for RNA that has folded into different regions.
 * @author Luke Czapla
 * @since 5.0.0-snapshot
 *
 */
public class TertiaryBasePairParameters extends BasePairParameters {

    // These are the criteria used to select proper base pairs.
    protected static double MaxStagger = 2.0, MaxPropeller = 60.0;

    public TertiaryBasePairParameters(Structure structure, boolean RNA, boolean removeDups) {
        super(structure, RNA, removeDups);
    }

    /**
     * This is an alternative implementation of findPair() that looks for anything that would fit the
     * criteria for a base-pair, useful for the context of tertiary structure of RNA.  Intra-strand base pairs
     * are found with this algorithm.
     * @param chains The list of chains already found to be nucleic acids
     * @return A list of the Pair of groups that match the base pair criteria, including intra-strand groups.
     */
    @Override
    public List<Pair<Group>> findPairs(List<Chain> chains) {
        List<Pair<Group>> result = new ArrayList<>();
        boolean lastFoundPair = false;
        for (int i = 0; i < chains.size(); i++) {
            Chain c = chains.get(i);
            String sequence = c.getAtomSequence();
            Integer type1, type2;
            for (int j = 0; j < sequence.length(); j++) {
                boolean foundPair = false;
                for (int k = sequence.length()-1; k >= j + 3 && !foundPair; k--) {
                    Group g1 = c.getAtomGroup(j);
                    Group g2 = c.getAtomGroup(k);
                    type1 = map.get(g1.getPDBName());
                    type2 = map.get(g2.getPDBName());
                    if (type1 == null || type2 == null) continue;
                    Atom a1 = g1.getAtom("C1'");
                    Atom a2 = g2.getAtom("C1'");
                    if (a1 == null || a2 == null) continue;
                    // C1'-C1' distance is one useful criteria
                    if (Math.abs(a1.getCoordsAsPoint3d().distance(a2.getCoordsAsPoint3d())-10.0) > 4.0) continue;
                    Pair<Group> ga = new Pair<>(g1, g2);
                    Matrix4d data = basePairReferenceFrame(ga);
                    // if the stagger is greater than 2 Å, it's not really paired.
                    if (Math.abs(pairParameters[5]) > MaxStagger) continue;
                    // if the propeller is ridiculous it's also not that good of a pair.
                    if (Math.abs(pairParameters[1]) > MaxPropeller) {
                        continue;
                    }
                    result.add(ga);
                    pairingNames.add(useRNA ? baseListRNA[type1]+baseListRNA[type2]: baseListDNA[type1]+baseListDNA[type2]);
                    foundPair = true;
                }
                if (!foundPair && lastFoundPair) {
                    if (pairSequence.length() > 0 && pairSequence.charAt(pairSequence.length()-1) != ' ')
                        pairSequence += ' ';
                }
                if (foundPair) pairSequence += (c.getAtomSequence().charAt(j));
                lastFoundPair = foundPair;
            }
        }
        result.addAll(super.findPairs(chains));
        return result;
    }

    public static double getMaxStagger() {
        return MaxStagger;
    }

    public static void setMaxStagger(double maxStagger) {
        MaxStagger = maxStagger;
    }

    public static double getMaxPropeller() {
        return MaxPropeller;
    }

    public static void setMaxPropeller(double maxPropeller) {
        MaxPropeller = maxPropeller;
    }
}
