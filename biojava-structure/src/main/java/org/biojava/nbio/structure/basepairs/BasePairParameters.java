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

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.geometry.SuperPosition;
import org.biojava.nbio.structure.geometry.SuperPositionQCP;
import org.biojava.nbio.structure.io.PDBFileReader;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.atan2;
import static java.lang.Math.acos;
import static java.lang.Math.PI;

/**
 * This module calculates the el Hassan-Calladine Base Pairing and Base-pair Step Parameters for any nucleic
 * acid containing structure that has the information about the core base-pair rings.
 * Citation: https://www.ncbi.nlm.nih.gov/pubmed/11601858
 *
 * The method that is usually overridden is findPairs(), this base implementation is used for a large-scale
 * analysis of the most proper helical regions in almost 4000 protein-DNA structures, almost
 * 2000 structures containing only DNA, or almost 1300 structures containing only RNA. (as of 7/2017).
 * Those who study tertiary structures for RNA folding should use the TertiaryBasePairParameters,
 * because this base class is only looking for base pairs between separate strands that exactly line up.
 * To relax the lining up policy and allow for non-canonical base pairs, use the MismatchedBasePairParameters
 * class, which will not consider intra-strand base pairing.
 *
 * @author Luke Czapla
 * @since 5.0.0
 *
 */
public class BasePairParameters implements Serializable {

    private static final long serialVersionUID = 6214502385L;
    private static Logger log = LoggerFactory.getLogger(BasePairParameters.class);

    // See URL http://ndbserver.rutgers.edu/ndbmodule/archives/reports/tsukuba/Table1.html
    // and the paper cited at the top of this class (also as Table 1).
    // These are hard-coded to avoid problems with resource paths.
    public static final String[] STANDARD_BASES = new String[] {
            "SEQRES   1 A    1  A\n" +
                    "ATOM      2  N9    A A   1      -1.291   4.498   0.000\n" +
                    "ATOM      3  C8    A A   1       0.024   4.897   0.000\n" +
                    "ATOM      4  N7    A A   1       0.877   3.902   0.000\n" +
                    "ATOM      5  C5    A A   1       0.071   2.771   0.000\n" +
                    "ATOM      6  C6    A A   1       0.369   1.398   0.000\n" +
                    "ATOM      8  N1    A A   1      -0.668   0.532   0.000\n" +
                    "ATOM      9  C2    A A   1      -1.912   1.023   0.000\n" +
                    "ATOM     10  N3    A A   1      -2.320   2.290   0.000\n" +
                    "ATOM     11  C4    A A   1      -1.267   3.124   0.000\n" +
                    "END",
            "SEQRES   1 A    1  G\n" +
                    "ATOM      2  N9    G A   1      -1.289   4.551   0.000\n" +
                    "ATOM      3  C8    G A   1       0.023   4.962   0.000\n" +
                    "ATOM      4  N7    G A   1       0.870   3.969   0.000\n" +
                    "ATOM      5  C5    G A   1       0.071   2.833   0.000\n" +
                    "ATOM      6  C6    G A   1       0.424   1.460   0.000\n" +
                    "ATOM      8  N1    G A   1      -0.700   0.641   0.000\n" +
                    "ATOM      9  C2    G A   1      -1.999   1.087   0.000\n" +
                    "ATOM     11  N3    G A   1      -2.342   2.364   0.001\n" +
                    "ATOM     12  C4    G A   1      -1.265   3.177   0.000\n" +
                    "END",
            "SEQRES   1 A    1  T\n" +
                    "ATOM      2  N1    T A   1      -1.284   4.500   0.000\n" +
                    "ATOM      3  C2    T A   1      -1.462   3.135   0.000\n" +
                    "ATOM      5  N3    T A   1      -0.298   2.407   0.000\n" +
                    "ATOM      6  C4    T A   1       0.994   2.897   0.000\n" +
                    "ATOM      8  C5    T A   1       1.106   4.338   0.000\n" +
                    "ATOM     10  C6    T A   1      -0.024   5.057   0.000\n" +
                    "END",
            "SEQRES   1 A    1  C\n" +
                    "ATOM      2  N1    C A   1      -1.285   4.542   0.000\n" +
                    "ATOM      3  C2    C A   1      -1.472   3.158   0.000\n" +
                    "ATOM      5  N3    C A   1      -0.391   2.344   0.000\n" +
                    "ATOM      6  C4    C A   1       0.837   2.868   0.000\n" +
                    "ATOM      8  C5    C A   1       1.056   4.275   0.000\n" +
                    "ATOM      9  C6    C A   1      -0.023   5.068   0.000\n" +
                    "END",
            "SEQRES   1 A    1  U\n" +
                    "ATOM      2  N1    U A   1      -1.284   4.500   0.000\n" +
                    "ATOM      3  C2    U A   1      -1.462   3.131   0.000\n" +
                    "ATOM      5  N3    U A   1      -0.302   2.397   0.000\n" +
                    "ATOM      6  C4    U A   1       0.989   2.884   0.000\n" +
                    "ATOM      8  C5    U A   1       1.089   4.311   0.000\n" +
                    "ATOM      9  C6    U A   1      -0.024   5.053   0.000\n"
    };

    // this is also hard-coded data about standard WC base pairs for both DNA and RNA
    protected static final String[] BASE_LIST_DNA = {"A", "G", "T", "C"};
    protected static final String[] BASE_LIST_RNA = {"A", "G", "U", "C"};
    protected static final Map<String, Integer> BASE_MAP;
   // private static List<String> RNAspecific = Arrays.asList("U", "URA"),
   //        DNAspecific = Arrays.asList("DC", "C", "CYT");
    protected static final Map<Integer, List<String>> RING_MAP;
    static {
        BASE_MAP = new HashMap<>();
        BASE_MAP.put("DA", 0); BASE_MAP.put("ADE", 0); BASE_MAP.put("A", 0);
        BASE_MAP.put("DG", 1); BASE_MAP.put("GUA", 1); BASE_MAP.put("G", 1);
        BASE_MAP.put("DT", 2); BASE_MAP.put("THY", 2); BASE_MAP.put("T", 2); BASE_MAP.put("U", 2); BASE_MAP.put("URA", 2);
        BASE_MAP.put("DC", 3); BASE_MAP.put("CYT", 3); BASE_MAP.put("C", 3);

        RING_MAP = new HashMap<>();
        RING_MAP.put(0, Arrays.asList("C8", "C2", "N3", "C4", "C5", "C6", "N7", "N1", "N9"));
        RING_MAP.put(1, Arrays.asList("C8", "C2", "N3", "C4", "C5", "C6", "N7", "N1", "N9"));
        RING_MAP.put(2, Arrays.asList("C6", "C2", "N3", "C4", "C5", "N1"));
        RING_MAP.put(3, Arrays.asList("C6", "C2", "N3", "C4", "C5", "N1"));
   }

    protected Structure structure;
    protected boolean canonical = true;
    protected boolean useRNA = false;
    protected boolean nonredundant = false;
    protected double[] pairParameters;

    // this is the main data that the user wants to get back out from the procedure.
    protected String pairSequence = "";
    protected double[][] pairingParameters;
    protected double[][] stepParameters;
    protected List<String> pairingNames = new ArrayList<>();
    protected List<Matrix4d> referenceFrames = new ArrayList<>();


    /**
     * This constructor takes a Structure object, finds base pair and base-pair step parameters
     * for double-helical regions within the structure.
     * @param structure The already-loaded structure to analyze.
     * @param useRNA whether to look for canonical RNA pairs.  By default (false) it analyzes DNA.
     * @param removeDups whether to only look for base-pair parameters for each unique sequence in
     *  the structure (if set to <i>true</i>)
     * @param canonical Whether to consider only Watson-Crick base pairs
     */
    public BasePairParameters(Structure structure, boolean useRNA, boolean removeDups, boolean canonical) {
        this.structure = structure;
        this.useRNA = useRNA;
        this.canonical = canonical;
        this.nonredundant = removeDups;

    }

    /**
     * This constructor takes a Structure object, whether to use RNA, and whether to remove duplicate sequences.
     * @param structure The already-loaded structure to analyze.
     * @param useRNA if true, the RNA standard bases will be used.  Otherwise, if false, it will work on standard DNA bases.
     * @param removeDups if true, duplicate sequences will not be considered.  This is for the analysis of X-ray structures from
     *                   RCSB, where there may be identical or similar units.
     */
    public BasePairParameters(Structure structure, boolean useRNA, boolean removeDups) {
        this(structure, useRNA, removeDups, false);
    }

    /**
     * This constructor takes a Structure object, and whether to use the RNA standard bases.
     * @param structure The already-loaded structure to analyze.
     * @param useRNA if true, the RNA standard bases will be used.  Otherwise, if false, it will work on standard DNA bases.
     */
    public BasePairParameters(Structure structure, boolean useRNA) {
        this(structure, useRNA, false, false);
    }

    /**
     * This constructor takes a Structure object, finds base pair and base-pair step parameters
     * for double-helical regions within the structure for only canonical DNA pairs.
     * @param structure The already-loaded structure to analyze.
     */
    public BasePairParameters(Structure structure) {
        this(structure, false, false, true);
    }


    /**
     * This method is the main function call to extract all step parameters, pairing parameters, and sequence
     * information from the Structure object provided to the constructor.
     * @return This same object with the populated data, convenient for output
     *  (e.g. <i>log.info(new BasePairParameters(structure).analyze());</i>)
     */
    public BasePairParameters analyze() {
        if (structure == null) {
            pairingParameters = null;
            stepParameters = null;
            return this;
        }
        List<Chain> nucleics = this.getNucleicChains(nonredundant);
        List<Pair<Group>> pairs = this.findPairs(nucleics);
        this.pairingParameters = new double[pairs.size()][6];
        this.stepParameters = new double[pairs.size()][6];
        Matrix4d lastStep;
        Matrix4d currentStep = null;
        for (int i = 0; i < pairs.size(); i++) {
            lastStep = currentStep;
            currentStep = this.basePairReferenceFrame(pairs.get(i));
            referenceFrames.add((Matrix4d)currentStep.clone());
            for (int j = 0; j < 6; j++) pairingParameters[i][j] = pairParameters[j];
            if (i != 0) {
                lastStep.invert();
                lastStep.mul(currentStep);
                double[] sparms = calculateTp(lastStep);
                for (int j = 0; j < 6; j++) stepParameters[i][j] = sparms[j];
            }
        }
        return this;
    }



    /**
     * This method returns the total number of base pairs that were found, used after the call to analyze().
     * @return An integer value, number of base pairs
     */
    public int getLength() {
        if (structure == null || pairParameters == null) throw new IllegalArgumentException("This structure is not analyzed or not initialized.");
        return pairingParameters.length;
    }


    /**
     * This method reports all the pair parameters, in the order of:
     * buckle, propeller, opening (in degrees), shear, stagger, stretch (in Å).
     * @return A double[][] with length equal to number of base pairs for rows, and 6 columns
     */
    public double[][] getPairingParameters() {
        return pairingParameters;
    }

    /**
     * This method reports all the base-pair step parameters, in the order of:
     * tilt, roll, twist (in degrees), shift, slide, rise (in Å).
     * @return A double[][] with length equal to number of base pairs (the first row 0 has no step
     *  and therefore is six zeroes), and 6 columns.
     */
    public double[][] getStepParameters() {
        return stepParameters;
    }


    /**
     * This method returns the primary strand's sequence where parameters were found.
     * There are spaces in the string anywhere there was a break in the helix or when
     * it goes from one helix to another helix in the structure. (the "step" is still returned)
     * @return String of primary sequence with spaces between gaps and new helices.
     */
    public String getPairSequence() {
        return pairSequence;
    }


    /**
     * This method returns the names of the pairs in terms of A, G, T/U, and C for each base pair group in the
     * list.  The first character is the leading strand base and the second character is the complementary base
     * @return
     */
    public List<String> getPairingNames() {
        return pairingNames;
    }

    public List<Matrix4d> getReferenceFrames() {
        return referenceFrames;
    }

    /**
     * This method is an internal test that the base pair specified is within a valid range.  If not, it throws an exception
     * with a message.
     * @param bp The index of the base pair or base-pair step to return.
     */
    private void checkArgument(int bp) {
        if (bp < 0 || bp >= getPairingParameters().length) throw new IllegalArgumentException("Base pair number is out of range.");
    }

    /**
     * This method returns the buckle in degrees for the given base pair
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in degrees)
     */
    public Double getBuckle(int bp) {
        checkArgument(bp);
        return pairingParameters[bp][0];
    }

    /**
     * This method returns the propeller ("propeller-twist") in degrees for the given base pair
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in degrees)
     */
    public Double getPropeller(int bp) {
        checkArgument(bp);
        return pairingParameters[bp][1];
    }

    /**
     * This method returns the opening in degrees for the given base pair
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in degrees)
     */
    public Double getOpening(int bp) {
        checkArgument(bp);
        return pairingParameters[bp][2];
    }

    /**
     * This method returns the shear in Å for the given base pair
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in Å)
     */
    public Double getShear(int bp) {
        checkArgument(bp);
        return pairingParameters[bp][3];
    }

    /**
     * This method returns the stretch in Å for the given base pair
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in Å)
     */
    public Double getStretch(int bp) {
        checkArgument(bp);
        return pairingParameters[bp][4];
    }

    /**
     * This method returns the stagger in Å for the given base pair
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in Å)
     */
    public Double getStagger(int bp) {
        checkArgument(bp);
        return pairingParameters[bp][5];
    }

    /**
     * This method returns the tilt for the given base pair, relative to the one before it.
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in degrees)
     */
    public Double getTilt(int bp) {
        checkArgument(bp);
        return stepParameters[bp][0];
    }

    /**
     * This method returns the roll for the given base pair, relative to the one before it.
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in degrees)
     */
    public Double getRoll(int bp) {
        if (bp < 0 || bp >= getStepParameters().length) throw new IllegalArgumentException("Base pair number is out of range.");
        return stepParameters[bp][1];
    }

    /**
     * This method returns the twist for the given base pair, relative to the one before it.
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in degrees)
     */
    public Double getTwist(int bp) {
        if (bp < 0 || bp >= getStepParameters().length) throw new IllegalArgumentException("Base pair number is out of range.");
        return stepParameters[bp][2];
    }

    /**
     * Return the shift for the given base pair, relative to the one before it.
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in Å)
     */
    public Double getShift(int bp) {
        if (bp < 0 || bp >= getStepParameters().length) throw new IllegalArgumentException("Base pair number is out of range.");
        return stepParameters[bp][3];
    }

    /**
     * This method returns the slide for the given base pair, relative to the one before it.
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in Å)
     */
    public Double getSlide(int bp) {
        if (bp < 0 || bp >= getStepParameters().length) throw new IllegalArgumentException("Base pair number is out of range.");
        return stepParameters[bp][4];
    }

    /**
     * This method returns the rise for the given base pair, relative to the one before it.
     * @param bp the number of the base pair (starting with 0)
     * @return the value as a double (in Å)
     */
    public Double getRise(int bp) {
        if (bp < 0 || bp >= getStepParameters().length) throw new IllegalArgumentException("Base pair number is out of range.");
        return stepParameters[bp][5];
    }


    /**
     * This method reports all the nucleic acid chains and has an option to remove duplicates if you
     * are considering an analysis of only unique DNA or RNA helices in the Structure.
     * @param removeDups If true, it will ignore duplicate chains
     * @return A list of all the nucleic acid chains in order of the Structure
     */
    public List<Chain> getNucleicChains(boolean removeDups) {
        if (structure == null) return new ArrayList<>();
        List<Chain> chains = structure.getChains();
        List<Chain> result = new ArrayList<>();
        for (Chain c: chains) {
            if (c.isNucleicAcid()) {
                result.add(c);
            }
        }
        if (removeDups) for (int i = 0; i < result.size(); i++) {
            for (int j = i+2; j < result.size(); j++) {
                // remove duplicate sequences (structures with two or more identical units)
                if (result.get(i).getAtomSequence().equals(result.get(j).getAtomSequence())) {
                    result.remove(j);
                }
            }
        }
        return result;
    }

    /**
     * This method performs a search for base pairs in the structure.  The criteria is alignment of
     * sequences and the canonical base pairs of DNA or RNA. Use MismatchedBasePairParameters
     * or TertiaryBasePairParameters for finding higher-order associations.
     * @param chains The list of chains already found to be nucleic acids
     * @return The list of corresponding Watson-Crick groups as pairs, as a Pair of nucleic acid Groups
     */
    public List<Pair<Group>> findPairs(List<Chain> chains) {
        List<Pair<Group>> result = new ArrayList<>();
        for (int i = 0; i < chains.size(); i++) {
            Chain c = chains.get(i);
            for (int j = i+1; j < chains.size(); j++) {
                String complement = complement(chains.get(j).getAtomSequence(), useRNA);
                String match = longestCommonSubstring(c.getAtomSequence(), complement);
                if (log.isDebugEnabled()) {
                    log.debug(c.getAtomSequence() + " " + chains.get(j).getAtomSequence() + " " + match);
                }
                int index1 = c.getAtomSequence().indexOf(match);
                int index2 = complement.length() - complement.indexOf(match) - 1;
                for (int k = 0; k < match.length(); k++) {
                    Group g1 = c.getAtomGroup(index1+k);
                    Group g2 = chains.get(j).getAtomGroup(index2-k);
                    Integer type1 = BASE_MAP.get(g1.getPDBName());
                    Integer type2 = BASE_MAP.get(g2.getPDBName());
                    if (type1 == null || type2 == null) {
                        if (pairSequence.length() != 0 && pairSequence.charAt(pairSequence.length()-1) != ' ') pairSequence += ' ';
                        continue;
                    }
                    Atom a1 = g1.getAtom(RING_MAP.get(type1).get(0));
                    Atom a2 = g2.getAtom(RING_MAP.get(type2).get(0));

                    if (a1 == null) {
                        log.info("Error processing " + g1.getPDBName());
                        if (pairSequence.length() != 0 && pairSequence.charAt(pairSequence.length()-1) != ' ') pairSequence += ' ';
                        continue;
                    }
                    if (a2 == null) {
                        log.info("Error processing " + g2.getPDBName());
                        if (pairSequence.length() != 0 && pairSequence.charAt(pairSequence.length()-1) != ' ') pairSequence += ' ';
                        continue;
                    }

                    double dx = a1.getX()-a2.getX();
                    double dy = a1.getY()-a2.getY();
                    double dz = a1.getZ()-a2.getZ();
                    double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);
                    //log.info("C8-C6 Distance (Å): " + distance);
                    // could be a base pair
                    if (Math.abs(distance-10.0) < 4.0) {
                        boolean valid = true;
                        for (String atomname : RING_MAP.get(type1)) {
                            Atom a = g1.getAtom(atomname);
                            if (a == null) valid = false;
                        }
                        if (valid) for (String atomname: RING_MAP.get(type2)) {
                            Atom a = g2.getAtom(atomname);
                            if (a == null) valid = false;
                        }
                        if (valid) {
                            result.add(new Pair<Group>(g1, g2));
                            pairingNames.add((useRNA ? BASE_LIST_RNA[type1]+ BASE_LIST_RNA[type2] : BASE_LIST_DNA[type1]+ BASE_LIST_DNA[type2]));
                            pairSequence += c.getAtomSequence().charAt(index1 + k);
                        } else if (pairSequence.length() != 0 && pairSequence.charAt(pairSequence.length()-1) != ' ') pairSequence += ' ';
                    } else if (pairSequence.length() != 0 && pairSequence.charAt(pairSequence.length()-1) != ' ') pairSequence += ' ';
                }
                if (pairSequence.length() != 0 && pairSequence.charAt(pairSequence.length()-1) != ' ') pairSequence += ' ';
            }
            //log.info();
        }
        log.info("Matched: " + pairSequence);
        return result;
    }


    /**
     * This method calculates the central frame (4x4 transformation matrix) of a single base pair.
     * @param pair An array of the two groups that make a hypothetical pair
     * @return The middle frame of the center of the base-pair formed
     */
    public Matrix4d basePairReferenceFrame(Pair<Group> pair) {
        Integer type1 = BASE_MAP.get(pair.getFirst().getPDBName());
        Integer type2 = BASE_MAP.get(pair.getSecond().getPDBName());
        SuperPosition sp = new SuperPositionQCP(true);
        if (type1 == null || type2 == null) return null;
        PDBFileReader pdbFileReader = new PDBFileReader();
        Structure s1, s2;
        try {
            s1 = pdbFileReader.getStructure(new ByteArrayInputStream(STANDARD_BASES[type1].getBytes()));
            s2 = pdbFileReader.getStructure(new ByteArrayInputStream(STANDARD_BASES[type2].getBytes()));
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        Group std1 = s1.getChain("A").getAtomGroup(0);
        Group std2 = s2.getChain("A").getAtomGroup(0);

        Point3d[] pointref = new Point3d[std1.getAtoms().size()];
        Point3d[] pointact = new Point3d[std1.getAtoms().size()];
        int count = 0;

        for (Atom a : std1.getAtoms()) {
            if (pair.getFirst().getAtom(a.getName()) == null) return null;
            pointref[count] = a.getCoordsAsPoint3d();
            pointact[count] = pair.getFirst().getAtom(a.getName()).getCoordsAsPoint3d();
            count++;
        }
        assert count == std1.getAtoms().size();
        Matrix4d ref1 = (Matrix4d)sp.superposeAndTransform(pointact, pointref).clone();

        pointref = new Point3d[std2.getAtoms().size()];
        pointact = new Point3d[std2.getAtoms().size()];

        count = 0;
        for (Atom a : std2.getAtoms()) {
            if (pair.getSecond().getAtom(a.getName()) == null) return null;
            pointref[count] = a.getCoordsAsPoint3d();
            pointact[count] = pair.getSecond().getAtom(a.getName()).getCoordsAsPoint3d();
            count++;
        }
        assert count == std2.getAtoms().size();

        Matrix4d temp = (Matrix4d)ref1.clone();
        Matrix4d temp2 = (Matrix4d)temp.clone();
        Matrix4d ref2 = sp.superposeAndTransform(pointact, pointref);

        double[][] v = new double[3][4];
        double[] y3 = new double[4];
        double[] z3 = new double[4];
        ref2.getColumn(1, y3);
        ref2.getColumn(2, z3);
        double[] z31 = new double[4];
        ref1.getColumn(2, z31);
        if (z3[0]*z31[0]+z3[1]*z31[1]+z3[2]*z31[2] < 0.0) {
            for (int i = 0; i < 3; i++) {
                y3[i] *= -1.0;
                z3[i] *= -1.0;
            }
        }
        ref2.setColumn(1, y3);
        ref2.setColumn(2, z3);

        temp.add(ref2);
        temp.mul(0.5);
        double[] x3 = new double[4];
        temp.getColumn(0, x3);
        temp.getColumn(1, y3);
        temp.getColumn(2, z3);
        x3 = removeComponent(x3, z3);
        x3 = removeComponent(x3, y3);
        y3 = removeComponent(y3, z3);
        temp.setColumn(0, x3);
        temp.setColumn(1, y3);
        temp.setColumn(2, z3);

        // normalize the short, long, and normal axes
        for (int i = 0; i < 3; i++) {
            temp.getColumn(i, v[i]);
            double r = Math.sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
            for (int j = 0; j < 3; j++) {
                v[i][j] /= r;
            }
            temp.setColumn(i, v[i]);
        }

        // calculate pairing parameters: buckle, propeller, opening, shear, stretch, stagger
        temp2.invert();
        temp2.mul(ref2);
        pairParameters = calculateTp(temp2);
        for (int i = 0; i < 6; i++) pairParameters[i] *= -1;

        // return the central frame of the base pair
        return temp;

    }


    @Override
    public String toString() {
        if (getPairingParameters() == null) return "No data";
        StringBuilder result = new StringBuilder(10000);
        result.append(pairingParameters.length + " base pairs\n");
        result.append("bp: buckle propeller opening shear stretch stagger tilt roll twist shift slide rise\n");
        for (int i = 0; i < pairingParameters.length; i++) {
            result.append(pairingNames.get(i)+": ");
            for (int j = 0; j < 6; j++)
                result.append(String.format("%5.4f", pairingParameters[i][j]) + " ");
            for (int j = 0; j < 6; j++)
                result.append(String.format("%5.4f", stepParameters[i][j]) + " ");
            result.append("\n");
        }
        return result.toString();
    }


    // The following methods are just helper classes for the rapid analyze of base-pair geometry.
    /**
     * This method calculates pairing and step parameters from 4x4 transformation matrices (used internally)
     * that comes out as a Matrix4d.
     * @param input the 4x4 matrix representing the transformation from strand II -> strand I or pair i to pair i+1
     * @return Six parameters as double[6]
     */
    public static double[] calculateTp(Matrix4d input) {

        double[][] A = new double[4][4];
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
            A[i][j] = input.getElement(i, j);
        }
        double[] M = new double[6];

        double cosgamma, gamma, phi, omega, sgcp, omega2_minus_phi,
                sm, cm, sp, cp, sg, cg;

        cosgamma = A[2][2];
        if (cosgamma > 1.0) cosgamma = 1.0;
        else if (cosgamma < -1.0) cosgamma = -1.0;

        gamma = acos(cosgamma);

        sgcp = A[1][1]*A[0][2]-A[0][1]*A[1][2];

        if (gamma == 0.0) omega = -atan2(A[0][1],A[1][1]);
        else omega = atan2(A[2][1]*A[0][2]+sgcp*A[1][2],sgcp*A[0][2]-A[2][1]*A[1][2]);

        omega2_minus_phi = atan2(A[1][2],A[0][2]);

        phi = omega/2.0 - omega2_minus_phi;

        M[0] = gamma*sin(phi)*180.0/PI;
        M[1] = gamma*cos(phi)*180.0/PI;
        M[2] = omega*180.0/PI;

        sm = sin(omega/2.0-phi);
        cm = cos(omega/2.0-phi);
        sp = sin(phi);
        cp = cos(phi);
        sg = sin(gamma/2.0);
        cg = cos(gamma/2.0);

        M[3] = (cm*cg*cp-sm*sp)*A[0][3]+(sm*cg*cp+cm*sp)*A[1][3]-sg*cp*A[2][3];
        M[4] = (-cm*cg*sp-sm*cp)*A[0][3]+(-sm*cg*sp+cm*cp)*A[1][3]+sg*sp*A[2][3];
        M[5] = (cm*sg)*A[0][3]+(sm*sg)*A[1][3]+cg*A[2][3];

        return M;

    }

    /**
     * This method returns the complement of a base. (used internally)
     * @param base The letter of the base
     * @param RNA Whether it is RNA (if false, it is DNA)
     * @return The character representing the complement of the base
     */
    protected static char complementBase(char base, boolean RNA) {
        if (base == 'A' && RNA) return 'U';
        if (base == 'A') return 'T';
        if (base == 'T' && !RNA) return 'A';
        if (base == 'U' && RNA) return 'A';
        if (base == 'C') return 'G';
        if (base == 'G') return 'C';
        return ' ';
    }

    /**
     * Simple helper method for quickly checking the complement of a sequence, see also DNASequence nad RNASequence classes
     * for more extensively useful functions not used in this narrow context of structural biology of base pairs.  (Used internally)
     */
    private static String complement(String sequence, boolean RNA) {
        String result = "";
        for (int i = sequence.length() - 1; i >= 0; i--) {
            result += complementBase(sequence.charAt(i), RNA);
        }
        return result;
    }

    /**
     * This does a 3D Vector cross product of two vectors as double arrays. (used internally)
     *
     * @param a An array of length 3 or 4 (4th component is ignored)
     * @param b An array of length 3 or 4 (4th component is ignored)
     * @return The cross product of the vectors (just the first three components
     */
    @SuppressWarnings("unused")
	private static double[] cross(double[] a, double[] b) {
        assert a.length >= 3 && b.length >= 3;
        double[] result = new double[4];
        result[0] = a[1]*b[2]-a[2]*b[1];
        result[1] = a[2]*b[0]-a[0]*b[2];
        result[2] = a[0]*b[1]-a[1]*b[0];
        return result;
    }

    /**
     * This method removes any component of vector a that is along vector b. (used internally)
     * @param a The array (vector) to remove component from
     * @param b The component array (vector) to remove from the first
     * @return The original array a with any component along b removed from it.
     */
    private static double[] removeComponent(double[] a, double[] b) {
        double dot = 0;
        double[] result = new double[4];
        for (int i = 0; i < 3; i++) {
            dot += a[i]*b[i];
        }
        for (int i = 0; i < 3; i++) {
            result[i] = a[i]-dot*b[i];
        }
        return result;

    }

    /**
     * This method finds the longest common substring between two strings. (used internally)
     * @param s1 The first string
     * @param s2 The second string
     * @return The substring itself
     */
    private static String longestCommonSubstring(String s1, String s2) {
        int start = 0;
        int max = 0;
        for (int i = 0; i < s1.length(); i++) {
            for (int j = 0; j < s2.length(); j++) {
                int x = 0;
                while (s1.charAt(i + x) == s2.charAt(j + x)) {
                    x++;
                    if (((i + x) >= s1.length()) || ((j + x) >= s2.length())) break;
                }
                if (x > max) {
                    max = x;
                    start = i;
                }
            }
        }
        return s1.substring(start, (start + max));
    }

    /**
     * This returns true if a is the complement of b, false otherwise. (used internally)
     * @param a First letter
     * @param b Potential matching letter
     * @param RNA Whether it is RNA (if false, DNA rules are used)
     * @return True if the bases are complementary.
     */
    protected static boolean match(char a, char b, boolean RNA) {
        if (a == 'A' && b == 'T' && !RNA) return true;
        if (a == 'A' && b == 'U' && RNA) return true;
        if (a == 'T' && b == 'A' && !RNA) return true;
        if (a == 'U' && b == 'A' && RNA) return true;
        if (a == 'G' && b == 'C') return true;
        if (a == 'C' && b == 'G') return true;
        return false;
    }
}
