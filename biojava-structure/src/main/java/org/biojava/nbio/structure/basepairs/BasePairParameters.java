package org.biojava.nbio.structure.basepairs;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.geometry.SuperPosition;
import org.biojava.nbio.structure.geometry.SuperPositionQCP;
import org.biojava.nbio.structure.io.PDBFileReader;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.*;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static java.lang.Math.*;

/**
 * Contributed to BioJava under its LGPL
 * This module calculates the el Hassan-Calladine Base Pairing and Base-pair Step Parameters
 * Citation: https://www.ncbi.nlm.nih.gov/pubmed/11601858
 *
 * The method that might be overridden is findPairs(), this implementation is used for a large-scale
 * analysis of the most proper helical regions in almost 4000 protein-DNA structures, almost
 * 2000 structures containing only DNA, or almost 1300 structures containing only RNA. (as of 7/2017).
 * Those who study tertiary structures for RNA folding would be better using their own method,
 * because this is only looking for base pairs between separate strands.
 *
 * Created by luke czapla on 7/20/17.
 */
public class BasePairParameters {

    private static Logger log = LoggerFactory.getLogger(BasePairParameters.class);

    // See URL http://ndbserver.rutgers.edu/ndbmodule/archives/reports/tsukuba/Table1.html
    // and the paper cited at the top of this class (also as Table 1).
    // These are hard-coded to avoid problems with resource paths.
    private static String[] standardBases = new String[] {
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
    protected static String[] baseListDNA = {"A", "G", "T", "C"};
    protected static String[] baseListRNA = {"A", "G", "U", "C"};
    protected static Map<String, Integer> map;
   // private static List<String> RNAspecific = Arrays.asList("U", "URA"),
   //        DNAspecific = Arrays.asList("DC", "C", "CYT");
    protected static Map<Integer, List<String>> ringMap;
    static {
        map = new HashMap<>();
        map.put("DA", 0); map.put("ADE", 0); map.put("A", 0);
        map.put("DG", 1); map.put("GUA", 1); map.put("G", 1);
        map.put("DT", 2); map.put("THY", 2); map.put("T", 2); map.put("U", 2); map.put("URA", 2);
        map.put("DC", 3); map.put("CYT", 3); map.put("C", 3);
        // chemically modified bases, leaving out to ignore (to treat as gaps) right now.
        //map.put("DZM", 0);
        //map.put("UCL", 2);
        //map.put("2DT", 2);
        //map.put("1CC", 3); map.put("5CM", 3);
        ringMap = new HashMap<>();
        ringMap.put(0, Arrays.asList("C8", "C2", "N3", "C4", "C5", "C6", "N7", "N1", "N9"));
        ringMap.put(1, Arrays.asList("C8", "C2", "N3", "C4", "C5", "C6", "N7", "N1", "N9"));
        ringMap.put(2, Arrays.asList("C6", "C2", "N3", "C4", "C5", "N1"));
        ringMap.put(3, Arrays.asList("C6", "C2", "N3", "C4", "C5", "N1"));
   }

    protected Structure structure;
    protected boolean useRNA = false;
    protected double[] pairParameters;

    // this is the main data that you want to get back out from the procedure.
    protected String pairSequence = "";
    protected double[][] pairingParameters;
    protected double[][] stepParameters;
    protected List<String> pairingNames = new ArrayList<>();


    /**
     * Constructor takes a Structure object, finds base pair and base-pair step parameters
     * for double-helical regions within the structure.
     * @param structure The already-loaded structure to analyze.
     * @param useRNA whether to look for canonical RNA pairs.  By default (false) it analyzes DNA.
     * @param removeDups whether to only look for base-pair parameters for each unique sequence in
     *  the structure (if set to <i>true</i>)
     */
    public BasePairParameters(Structure structure, boolean useRNA, boolean removeDups) {
        this.structure = structure;
        this.useRNA = useRNA;
        if (structure == null) {
            pairingParameters = null;
            stepParameters = null;
            return;
        }
        List<Chain> nucleics = this.getNucleicChains(removeDups);
        List<Group[]> pairs = this.findPairs(nucleics);
        pairingParameters = new double[pairs.size()][6];
        stepParameters = new double[pairs.size()][6];
        Matrix4d lastStep = null;
        Matrix4d currentStep = null;
        for (int i = 0; i < pairs.size(); i++) {
            lastStep = currentStep;
            currentStep = this.basePairReferenceFrame(pairs.get(i));
            for (int j = 0; j < 6; j++) pairingParameters[i][j] = pairParameters[j];
            if (i != 0) {
                lastStep.invert();
                lastStep.mul(currentStep);
                double[] sparms = calculatetp(lastStep);
                for (int j = 0; j < 6; j++) stepParameters[i][j] = sparms[j];
            }
;       }

    }


    /**
     * Constructor takes a Structure object, finds base pair and base-pair step parameters
     * for double-helical regions within the structure for only canonical DNA pairs.
     * @param structure The already-loaded structure to analyze.
     */
    public BasePairParameters(Structure structure) {
        this(structure, false, false);
    }

    /**
     * This reports all the pair parameters, in the order of:
     * buckle, propeller, opening (in degrees), shear, stagger, stretch (in Å).
     * @return A double[][] with length equal to number of base pairs for rows, and 6 columns
     */
    public double[][] getPairingParameters() {
        return pairingParameters;
    }

    /**
     * This reports all the base-pair step parameters, in the order of:
     * tilt, roll, twist (in degrees), shift, slide, rise (in Å).
     * @return A double[][] with length equal to number of base pairs (the first row 0 has no step
     *  and therefore is six zeroes), and 6 columns.
     */
    public double[][] getStepParameters() {
        return stepParameters;
    }

    /**
     * This returns the primary strand's sequence where parameters were found.
     * There are spaces in the string anywhere there was a break in the helix or when
     * it goes from one helix to another helix in the structure. (the "step" is still returned!)
     * @return String of primary sequence with spaces between gaps and new helices.
     */
    public String getPairSequence() {
        return pairSequence;
    }

    /**
     * This returns the names of the pairs in terms of A, G, T/U, and C for each base pair group in the
     * list.  The first character is the leading strand base and the second character is the complementary base
     * @return
     */
    public List<String> getPairingNames() {
        return pairingNames;
    }

    /**
     * This reports all the nucleic acid chains and has an option to remove duplicates if you
     * are considering an analyze of only unique DNA or RNA helices in the Structure.
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
                if (result.get(i).getSeqResSequence().equals(result.get(j).getSeqResSequence())) {
                    result.remove(j);
                }
            }
        }
        return result;
    }


    /**
     * This performs a search for base pairs in the structure.  The criteria is alignment of
     * sequences and the canonical base pairs of DNA and RNA.
     * @param chains The list of chains already found to be nucleic acids
     * @return The list of corresponding Watson-Crick groups as pairs, element 0 is on the
     *  forward strand and element 1 is on the reverse strand.
     */
    public List<Group[]> findPairs(List<Chain> chains) {
        List<Group[]> result = new ArrayList<>();
        for (int i = 0; i < chains.size(); i++) {
            Chain c = chains.get(i);
            for (int j = i+1; j < chains.size(); j++) {
                String complement = complement(chains.get(j).getSeqResSequence(), useRNA);
                String match = longestCommonSubstring(c.getSeqResSequence(), complement);
                //log.info(c.getSeqResSequence() + " " + chains.get(j).getSeqResSequence() + " " + match);
                int index1 = c.getSeqResSequence().indexOf(match);
                int index2 = complement.length() - complement.indexOf(match) - 1;
                for (int k = 0; k < match.length(); k++) {
                    Group g1 = c.getSeqResGroup(index1+k);
                    Group g2 = chains.get(j).getSeqResGroup(index2-k);
                    Integer type1 = map.get(g1.getPDBName());
                    Integer type2 = map.get(g2.getPDBName());
                    if (type1 == null || type2 == null) {
                        if (pairSequence.length() != 0 && pairSequence.charAt(pairSequence.length()-1) != ' ') pairSequence += ' ';
                        continue;
                    }
                    Atom a1 = g1.getAtom(ringMap.get(type1).get(0));
                    Atom a2 = g2.getAtom(ringMap.get(type2).get(0));

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
                    if (Math.abs(distance-10.0) < 5.0) {
                        boolean valid = true;
                        for (String atomname : ringMap.get(type1)) {
                            Atom a = g1.getAtom(atomname);
                            if (a == null) valid = false;
                        }
                        if (valid) for (String atomname: ringMap.get(type2)) {
                            Atom a = g2.getAtom(atomname);
                            if (a == null) valid = false;
                        }
                        if (valid) {
                            result.add(new Group[]{g1, g2});
                            pairingNames.add((useRNA ? baseListRNA[type1]+baseListRNA[type2] : baseListDNA[type1]+baseListDNA[type2]));
                            pairSequence += c.getSeqResSequence().charAt(index1 + k);
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


    public Matrix4d basePairReferenceFrame(Group[] pair) {
        Integer type1 = map.get(pair[0].getPDBName());
        Integer type2 = map.get(pair[1].getPDBName());
        SuperPosition sp = new SuperPositionQCP(true);
        if (type1 == null || type2 == null) return null;
        PDBFileReader pdbFileReader = new PDBFileReader();
        Structure s1, s2;
        try {
            s1 = pdbFileReader.getStructure(new ByteArrayInputStream(standardBases[type1].getBytes()));
            s2 = pdbFileReader.getStructure(new ByteArrayInputStream(standardBases[type2].getBytes()));
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
            if (pair[0].getAtom(a.getName()) == null) return null;
            pointref[count] = a.getCoordsAsPoint3d();
            pointact[count] = pair[0].getAtom(a.getName()).getCoordsAsPoint3d();
            count++;
        }
        assert count == std1.getAtoms().size();
        Matrix4d ref1 = (Matrix4d)sp.superposeAndTransform(pointact, pointref).clone();

        pointref = new Point3d[std2.getAtoms().size()];
        pointact = new Point3d[std2.getAtoms().size()];

        count = 0;
        for (Atom a : std2.getAtoms()) {
            if (pair[1].getAtom(a.getName()) == null) return null;
            pointref[count] = a.getCoordsAsPoint3d();
            pointact[count] = pair[1].getAtom(a.getName()).getCoordsAsPoint3d();
            count++;
        }
        assert count == std2.getAtoms().size();

        //    log.info(ref1);
        Matrix4d temp = (Matrix4d)ref1.clone();
        Matrix4d temp2 = (Matrix4d)temp.clone();
        Matrix4d ref2 = sp.superposeAndTransform(pointact, pointref);
        //    log.info(ref2);
        double[][] v = new double[3][4];
        double[] y3 = new double[4];
        double[] z3 = new double[4];
        ref2.getColumn(1, y3);
        ref2.getColumn(2, z3);
        for (int i = 0; i < 3; i++) {
            y3[i] *= -1.0;
            z3[i] *= -1.0;
        }
        ref2.setColumn(1, y3);
        ref2.setColumn(2, z3);
        temp.add(ref2);
        temp.mul(0.5);
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
        pairParameters = calculatetp(temp2);
        for (int i = 0; i < 6; i++) pairParameters[i] *= -1;

        // return the central frame of the base pair
        return temp;

    }



    /**
     * This method calculates pairing and step parameters from 4x4 transformation matrices
     * that come out as Matrix4d;
     * @param input the 4x4 matrix representing the transformation from strand II -> strand I or pair i to pair i+1
     * @return Six parameters as double[6]
     */
    public static double[] calculatetp(Matrix4d input) {

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


    public static char complementBase(char base, boolean RNA) {
        if (base == 'A' && RNA) return 'U';
        if (base == 'A') return 'T';
        if (base == 'T' && !RNA) return 'A';
        if (base == 'U' && RNA) return 'A';
        if (base == 'C') return 'G';
        if (base == 'G') return 'C';
        return ' ';
    }

    public static String complement(String sequence, boolean RNA) {
        String result = "";
        for (int i = sequence.length() - 1; i >= 0; i--) {
            result += complementBase(sequence.charAt(i), RNA);
        }
        return result;
    }


    public static String longestCommonSubstring(String s1, String s2) {
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

}
