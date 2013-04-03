/*
 *                  BioJava development code
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
 * Created on May 21, 2006
 *
 */
package org.biojava.bio.structure.align;

import org.biojava.bio.structure.StructureTools;

/** A class that contains all the parameters of the structure alignment algorithm.
 *
 * @author Andreas Prlic
 * @since 1.5
 * @version %I% %G%
 */
public class StrucAligParameters {


    int initialK;
    String[] usedAtomNames = { StructureTools.caAtomName, } ;

    // step 1

    int seedFragmentLength; // seed fragment length
    float seedRmsdCutoff;

    int fragmentLength;
    int diagonalDistance;
    int diagonalDistance2; // set to < 1 if not used.
    boolean reduceInitialFragments;

    // step 2
    float fragmentMiniDistance;
    float fragCompat;   // fragment compatibility cutoff
    int maxrefine;      // max number of JointFragments to be refined

    boolean joinPlo; // joining according to BioPython variant
    boolean joinFast; // apply a fast procedure for extending the alignments

    // joininf of fragments - checks
    boolean doAngleCheck  ;
    boolean doDistanceCheck;
    boolean doDensityCheck;
    boolean doRMSCheck;
    float densityCutoff;

    int angleDiff;      // directional difference
    double joinRMSCutoff; // rms cutoff to be applied during joining of fragments.

    // step 4
    float create_co; //  alignment generation cutoff
    int maxIter; //  # max number of iterations in refinement
    float gapOpen;// gap open penalty
    float gapExtension; // gap extensionpenalty
    int permutationSize; // minimal size for a permutated fragment ( -1 means no circular permutation search)
    float evalCutoff; //  alignment evaluation cutoff

    public StrucAligParameters() {
        super();
        setDefault();
    }

    public static StrucAligParameters getDefaultParameters(){
        StrucAligParameters n = new StrucAligParameters();
        return n;
    }

    private void setDefault() {
        initialK = 6;

        // step 1
        seedFragmentLength      = 8;
        seedRmsdCutoff          = 3.0f; // orig 2.0 - better?
        fragmentLength          = 10;
        diagonalDistance        = 3;
        diagonalDistance2       = 9; // this helps a lot in 1buz vs 1aua
        fragmentMiniDistance    = 3.5f; // orig 2
        angleDiff               = 10;
        fragCompat              = 6.0f; // orig 4.0
        maxrefine               = 20; // orig 20

        // step 2
        reduceInitialFragments  =  true; // if this is disabled, you might want to also disable doRMSCheck for large structures...
        joinRMSCutoff           = 5.0; // orig 4
        joinPlo                 = false;
        joinFast                = false;


        // 3 joint fragments
        doAngleCheck            = true;
        doDistanceCheck         = true;
        doRMSCheck              = true;

        doDensityCheck          = false; // hm this one needs improvements before being used
        densityCutoff           = 7.0f;

        //  step 3
        create_co           = 6.0f;
        maxIter             = 4;  // number of times dynamic programming is run. set to zero for quick search (but imprecise)
        gapOpen             = 20.0f;
        gapExtension        = 0.0f;
        permutationSize     = 20;
        evalCutoff          = 6.0f;
    }
    public String toString() {
        StringBuffer buf = new StringBuffer();
        String t = " ";

        Object[] params = new Object[]{new Integer(initialK) ,new Integer(seedFragmentLength),
        		new Float(seedRmsdCutoff),
        		new Integer(fragmentLength),
                new Integer(diagonalDistance), new Integer(diagonalDistance2), new Float(fragmentMiniDistance),
                new Integer(angleDiff),
                new Float(fragCompat), new Integer(maxrefine),
                new Boolean(reduceInitialFragments), new Double(joinRMSCutoff), new Boolean(joinPlo),
                new Boolean(doAngleCheck), new Boolean(doDistanceCheck), new Boolean(doRMSCheck),
                new Boolean(doDensityCheck), new Float(densityCutoff), new Float(create_co), new Integer(maxIter),
                new Float(gapOpen), new Float(gapExtension), new Integer(permutationSize), new Float(evalCutoff)};

        for (int i=0 ; i< params.length ; i++){
            buf.append(params[i]);
            buf.append(t);
        }


        return buf.toString();

    }

    public static StrucAligParameters getDBSearchParameters(){
        StrucAligParameters params = new StrucAligParameters();

        params.setMaxIter(0); // not so nice alignments, but significant similarities should already be found,
        // one could do a second interation later over the top ranking hits and make a nicer alignment.


        return params;
    }

    public float getDensityCutoff() {
        return densityCutoff;
    }

    public void setDensityCutoff(float densityCutoff) {
        this.densityCutoff = densityCutoff;
    }

    public int getInitialK() {
        return initialK;
    }

    public void setInitialK(int initialK) {
        this.initialK = initialK;
    }

    public int getSeedFragmentLength() {
        return seedFragmentLength;
    }

    public boolean isJoinFast(){
       return joinFast;
    }

    public void setJoinFast(boolean fastJoin){
       joinFast = fastJoin;
    }

    public boolean isJoinPlo() {
        return joinPlo;
    }

    public void setJoinPlo(boolean joinPlo) {
        this.joinPlo = joinPlo;
    }

    public void setSeedFragmentLength(int seedFragmentLength) {
        this.seedFragmentLength = seedFragmentLength;
    }



    public float getSeedRmsdCutoff() {
        return seedRmsdCutoff;
    }



    public void setSeedRmsdCutoff(float seedRmsdCutoff) {
        this.seedRmsdCutoff = seedRmsdCutoff;
    }



    public boolean isDoAngleCheck() {
        return doAngleCheck;
    }

    public void setDoAngleCheck(boolean doAngleCheck) {
        this.doAngleCheck = doAngleCheck;
    }

    public boolean isDoDensityCheck() {
        return doDensityCheck;
    }

    public void setDoDensityCheck(boolean doDensityCheck) {
        this.doDensityCheck = doDensityCheck;
    }

    public boolean isDoDistanceCheck() {
        return doDistanceCheck;
    }

    public void setDoDistanceCheck(boolean doDistanceCheck) {
        this.doDistanceCheck = doDistanceCheck;
    }

    public boolean isDoRMSCheck() {
        return doRMSCheck;
    }

    public void setDoRMSCheck(boolean doRMSCheck) {
        this.doRMSCheck = doRMSCheck;
    }

    public double getJoinRMSCutoff() {
        return joinRMSCutoff;
    }



    public void setJoinRMSCutoff(double joinRMSCutoff) {
        this.joinRMSCutoff = joinRMSCutoff;
    }



    public float getEvalCutoff() {
        return evalCutoff;
    }

    public void setEvalCutoff(float evalCutoff) {
        this.evalCutoff = evalCutoff;
    }




    public int getPermutationSize() {
        return permutationSize;
    }

    public void setPermutationSize(int permutationSize) {
        this.permutationSize = permutationSize;
    }

    public float getGapExtension() {
        return gapExtension;
    }

    public void setGapExtension(float gapExtension) {
        this.gapExtension = gapExtension;
    }

    public float getGapOpen() {
        return gapOpen;
    }

    public void setGapOpen(float gapOpen) {
        this.gapOpen = gapOpen;
    }

    public int getMaxIter() {
        return maxIter;
    }

    public void setMaxIter(int maxIter) {
        this.maxIter = maxIter;
    }

    public float getCreate_co() {
        return create_co;
    }

    public void setCreate_co(float create_co) {
        this.create_co = create_co;
    }

    /** if this is set to false, the time spent to joint the initial fragments (step 2)
     * is increased. - particular for large structures this increases calc. time a lot.
     * advantage: more combinations of fragments are used.
     *
     * @return a flag if the inital fragments should be reduced
     */
    public boolean reduceInitialFragments() {
        return reduceInitialFragments;
    }

    public void setReduceInitialFragments(boolean reduceInitialFragments) {
        this.reduceInitialFragments = reduceInitialFragments;
    }

    public int getAngleDiff() {
        return angleDiff;
    }

    public void setAngleDiff(int angleDiff) {
        this.angleDiff = angleDiff;
    }

    public float getFragCompat() {
        return fragCompat;
    }

    public void setFragCompat(float fragCompat) {
        this.fragCompat = fragCompat;
    }

    public int getMaxrefine() {
        return maxrefine;
    }

    public void setMaxrefine(int maxrefine) {
        this.maxrefine = maxrefine;
    }

    public String[] getUsedAtomNames() {
        return usedAtomNames;
    }

    public void setUsedAtomNames(String[] usedAtomNames) {
        this.usedAtomNames = usedAtomNames;
    }

    public int getFragmentLength() {
        return fragmentLength;
    }

    public void setFragmentLength(int fragmentLength) {
        this.fragmentLength = fragmentLength;
    }

    public int getDiagonalDistance() {
        return diagonalDistance;
    }

    public void setDiagonalDistance(int diagonalDistance) {
        this.diagonalDistance = diagonalDistance;
    }



    public int getDiagonalDistance2() {
        return diagonalDistance2;
    }

    public void setDiagonalDistance2(int diagonalDistance2) {
        this.diagonalDistance2 = diagonalDistance2;
    }

    public float getFragmentMiniDistance() {
        return fragmentMiniDistance;
    }

    public void setFragmentMiniDistance(float fragmentMiniDistance) {
        this.fragmentMiniDistance = fragmentMiniDistance;
    }




}
