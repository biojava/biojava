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
 * This class allows for finding inter-strand base pairs that are not necessarily canonical Watson-Crick pairs.
 * The implementation of findPair is different than that of the base class.  This class does not consider intra-strand
 * base pairing and for that, the TertiaryBasePairParameters class should be used.
 * @author Luke Czapla
 * @since 5.0.0
 *
 */
public class MismatchedBasePairParameters extends BasePairParameters {

	private static final long serialVersionUID = 2837124340169886674L;
	
	public static final double DEFAULT_MAX_STAGGER = 2.0;
    public static final double DEFAULT_MAX_PROPELLER = 60.0;
    public static final double DEFAULT_MAX_SHEAR = 5.0;
    public static final double DEFAULT_MAX_STRETCH = 5.0;

    // These are the criteria used to select proper base pairs.
    private double maxStagger = DEFAULT_MAX_STAGGER,
            maxShear = DEFAULT_MAX_SHEAR,
            maxStretch = DEFAULT_MAX_STRETCH,
            maxPropeller = DEFAULT_MAX_PROPELLER;

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
                        // TODO is this call needed?? JD 2018-03-07
                        @SuppressWarnings("unused")
						Matrix4d data = basePairReferenceFrame(ga);
                        // if the stagger is greater than 2 Å, it's not really paired.
                        if (Math.abs(pairParameters[5]) > maxStagger) continue;
                        // similarly, extreme shear and stretch is not a good base pair
                        if (Math.abs(pairParameters[3]) > maxShear) continue;
                        if (Math.abs(pairParameters[4]) > maxStretch) continue;

                        // if the propeller is ridiculous it's also not that good of a pair.
                        if (Math.abs(pairParameters[1]) > maxPropeller) {
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

    /**
     * This method returns the maximum stagger between bases used as criteria for the characterization of two bases as being paired.
     * @return the maximum propeller ("propeller-twist", in degrees) allowed.
     */
    public double getMaxStagger() {
        return maxStagger;
    }

    /**
     * This method sets the maximum stagger allowed for a base pair, prior to analyze() call
     * @param maxStagger The maximum propeller (in Å) allowed to consider two bases paired
     */
    public void setMaxStagger(double maxStagger) {
        this.maxStagger = maxStagger;
    }

    /**
     * This method returns the maximum shear between bases used as criteria for the characterization of two bases as being paired.
     * @return the maximum shear (in Å) allowed.
     */
    public double getMaxShear() {
        return maxShear;
    }

    /**
     * This method sets the maximum shear allowed for a base pair, prior to analyze() call
     * @param maxShear The maximum shear (in Å) allowed to consider two bases paired
     */
    public void setMaxShear(double maxShear) {
        this.maxShear = maxShear;
    }

    /**
     * This method returns the maximum stretch between bases used as criteria for the characterization of two bases as being paired.
     * @return the maximum stretch (in Å) allowed.
     */
    public double getMaxStretch() {
        return maxStretch;
    }

    /**
     * This method sets the maximum stretch allowed for a base pair, prior to analyze() call.
     * @param maxStretch The maximum stretch (in Å) allowed to consider two bases paired
     */
    public void setMaxStretch(double maxStretch) {
        this.maxStretch = maxStretch;
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
