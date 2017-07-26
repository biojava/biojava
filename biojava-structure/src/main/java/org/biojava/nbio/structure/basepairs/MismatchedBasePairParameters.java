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
 * Contributed to BioJava under its LGPL
 * This class allows for finding inter-strand base pairs that are not necessarily canonical Watson-Crick pairs.
 * The implementation of findPair is different than that of the base class.
 * @author Luke Czapla
 * @since 5.0.0-snapshot
 *
 */
public class MismatchedBasePairParameters extends BasePairParameters {

    // These are the criteria used to select proper base pairs.
    protected static double MaxStagger = 2.0, MaxShear = 5.0, MaxStretch = 5.0,
            MaxPropeller = 60.0;

    public MismatchedBasePairParameters(Structure structure, boolean RNA, boolean removeDups, boolean canonical) {

        super(structure, RNA, removeDups, canonical);

    }

    /**
     * This is an implementation for finding non-canonical base pairs when there may be missing or overhanging bases.
     * @param chains The list of chains already found to be nucleic acids
     * @return The list of the atom groups (residues) that are pairs, as a Pair of nucleic acid Groups
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
                        type1 = map.get(g1.getPDBName());
                        type2 = map.get(g2.getPDBName());
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
                        if (Math.abs(pairParameters[3]) > MaxShear) continue;
                        if (Math.abs(pairParameters[4]) > MaxStretch) continue;

                        // if the propeller is ridiculous it's also not that good of a pair.
                        if (Math.abs(pairParameters[1]) > MaxPropeller) {
                            continue;
                        }
                        result.add(ga);
                        pairingNames.add(useRNA ? baseListRNA[type1] + baseListRNA[type2] : baseListDNA[type1] + baseListDNA[type2]);
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


}
