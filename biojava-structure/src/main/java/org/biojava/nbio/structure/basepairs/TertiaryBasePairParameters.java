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
import java.util.ArrayList;
import java.util.List;

/**
 * This class also finds the base pairing and base-pair step parameters but has a broader definition
 * of a base pair so that non-canonical-WC base pairs will be detected and reported.  This is useful
 * for RNA that has folded into different regions, and for higher-order DNA structures.  Intra-strand
 * pairings are considered in this class (but in not the base class or MismatchedBasePairParameters class)
 * @author Luke Czapla
 * @since 5.0.0
 *
 */
public class TertiaryBasePairParameters extends BasePairParameters {

	private static final long serialVersionUID = 2556427111533466577L;
	
	public static final double DEFAULT_MAX_STAGGER = 2.0;
    public static final double DEFAULT_MAX_PROPELLER = 60.0;
    // These are the criteria used to select proper base pairs.
    private double maxStagger = DEFAULT_MAX_STAGGER,
            maxPropeller = DEFAULT_MAX_PROPELLER;

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
                    type1 = BASE_MAP.get(g1.getPDBName());
                    type2 = BASE_MAP.get(g2.getPDBName());
                    if (type1 == null || type2 == null) continue;
                    Atom a1 = g1.getAtom("C1'");
                    Atom a2 = g2.getAtom("C1'");
                    if (a1 == null || a2 == null) continue;
                    // C1'-C1' distance is one useful criteria
                    if (Math.abs(a1.getCoordsAsPoint3d().distance(a2.getCoordsAsPoint3d())-10.0) > 4.0) continue;
                    Pair<Group> ga = new Pair<>(g1, g2);
                    // TODO is this call needed?? JD 2018-03-07
                    @SuppressWarnings("unused")
					Matrix4d data = basePairReferenceFrame(ga);
                    // if the stagger is greater than 2 Å, it's not really paired.
                    if (Math.abs(pairParameters[5]) > maxStagger) continue;
                    // if the propeller is ridiculous it's also not that good of a pair.
                    if (Math.abs(pairParameters[1]) > maxPropeller) {
                        continue;
                    }
                    result.add(ga);
                    pairingNames.add(useRNA ? BASE_LIST_RNA[type1]+ BASE_LIST_RNA[type2]: BASE_LIST_DNA[type1]+ BASE_LIST_DNA[type2]);
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

    /**
     * This method returns the maximum stagger between bases used as criteria for the characterization of two bases as being paired.
     * @return the maximum stagger (in Å) allowed.
     */
    public double getMaxStagger() {
        return maxStagger;
    }

    /**
     * This method sets the maximum stagger allowed for a base pair, prior to analyze() call
     * @param maxStagger The maximum stagger (in Å) allowed to consider two bases paired
     */
    public void setMaxStagger(double maxStagger) {
        this.maxStagger = maxStagger;
    }

    /**
     * This method returns the maximum propeller twist between bases used as criteria for the characterization of two bases as being paired.
     * @return the maximum propeller ("propeller-twist", in degrees) allowed.
     */
    public double getMaxPropeller() {
        return maxPropeller;
    }

    /**
     * This method sets the maximum propeller allowed for a base pair, prior to analyze() call
     * @param maxPropeller The maximum propeller ("propeller-twist", in degrees) allowed to consider two bases paired
     */
    public void setMaxPropeller(double maxPropeller) {
        this.maxPropeller = maxPropeller;
    }
}
