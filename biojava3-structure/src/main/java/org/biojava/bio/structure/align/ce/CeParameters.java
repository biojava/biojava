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
 * Created on Sep 15, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;

import java.util.ArrayList;
import java.util.List;

import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;


/** Contains the parameters that can be sent to CE
 * 
 * @author Andreas Prlic
 *
 */
public class CeParameters implements ConfigStrucAligParams  {

	protected int winSize;
	protected double rmsdThr;
	protected double rmsdThrJoin;
	protected double maxOptRMSD;

	public static final int DEFAULT_SCORING_STRATEGY = 0;
	public static final int SIDE_CHAIN_SCORING = 1;
	public static final int SIDE_CHAIN_ANGLE_SCORING = 2;
	public static final int CA_AND_SIDE_CHAIN_ANGLE_SCORING = 3;
	public static final int SEQUENCE_CONSERVATION = 4;

	public static final String SCORING_STRATEGY = "ScoringStrategy";
	protected int scoringStrategy;
	//String[] alignmentAtoms;
	protected int maxGapSize;

	protected boolean showAFPRanges;
	protected int  sideChainScoringType;

	protected static final double DEFAULT_GAP_OPEN = 5.0;
	protected static final double DEFAULT_GAP_EXTENSION  = 0.5;
	protected static final double DISTANCE_INCREMENT=0.5;
	protected static final double DEFAULT_oRmsdThr = 2.0; 

	protected double gapOpen;
	protected double gapExtension;
	protected double distanceIncrement;
	protected double oRmsdThr;

	protected int maxNrIterationsForOptimization;

	protected SubstitutionMatrix<AminoAcidCompound> substitutionMatrix;	
	protected double seqWeight;

	public CeParameters(){
		reset();
	}



	@Override
	public String toString() {
		return "CeParameters [scoringStrategy=" + scoringStrategy 
		+ ", maxGapSize=" + maxGapSize 
		+ ", rmsdThr=" + rmsdThr 
		+ ", rmsdThrJoin="+ rmsdThrJoin 
		+ ", winSize=" + winSize 
		+ ", showAFPRanges=" + showAFPRanges 
		+ ", maxOptRMSD=" + maxOptRMSD
		+ ", seqWeight=" + seqWeight
		+ "]";
	}



	public void reset(){
		winSize = 8;
		rmsdThr = 3.0;
		rmsdThrJoin = 4.0;
		scoringStrategy = DEFAULT_SCORING_STRATEGY;
		maxGapSize = 30;
		showAFPRanges = false;
		maxOptRMSD = 99;

		gapOpen = DEFAULT_GAP_OPEN;
		gapExtension = DEFAULT_GAP_EXTENSION;
		distanceIncrement = DISTANCE_INCREMENT;
		oRmsdThr = DEFAULT_oRmsdThr;

		maxNrIterationsForOptimization = Integer.MAX_VALUE;
		seqWeight = 0;
	}

	/** The window size to look at
	 * 
	 * @return window size
	 */
	public Integer getWinSize() {
		return winSize;
	}
	public void setWinSize(Integer winSize) {
		this.winSize = winSize;
	}

	/** RMSD Threshold
	 * 
	 * @return RMSD threshold
	 */
	public Double getRmsdThr() {
		return rmsdThr;
	}
	public void setRmsdThr(Double rmsdThr) {
		this.rmsdThr = rmsdThr;
	}

	/** RMSD threshold for joining of AFPs
	 * 
	 * @return rmsd threshold
	 */
	public Double getRmsdThrJoin() {
		return rmsdThrJoin;
	}
	public void setRmsdThrJoin(Double rmsdThrJoin) {
		this.rmsdThrJoin = rmsdThrJoin;
	}

	public Integer getScoringStrategy()
	{
		return scoringStrategy;
	}


	/** Set the scoring strategy to use. 0 is default CE scoring scheme. 1 uses
	 * Side chain orientation.
	 * 
	 * @param scoringStrategy
	 */
	public void setScoringStrategy(Integer scoringStrategy)
	{
		this.scoringStrategy = scoringStrategy;
	}



	/** Set the Max gap size parameter. Default 30. For unlimited gaps set to -1
	 * 
	 * @param maxGapSize
	 */
	public void setMaxGapSize(Integer maxGapSize){
		this.maxGapSize = maxGapSize;
	}

	/** the Max gap size parameter G . default is 30, which was
	 * described to obtained empirically in the CE paper.
	 * the larger the max gap size, the longer the compute time,
	 * but in same cases drastically improved results. Set to -1 for unlimited gap size. 
	 * 
	 * @return max gap size parameter
	 */
	public Integer getMaxGapSize() {
		return maxGapSize;
	}


	public List<String> getUserConfigHelp() {
		List<String> params =new ArrayList<String>();
		String helpMaxGap = "This parameter configures the maximum gap size G, that is applied during the AFP extension. The larger the value, the longer the calculation time can become, Default value is 30. Set to 0 for no limit. " ;
		//String helpRmsdThr = "This configures the RMSD threshold applied during the trace of the fragment matrix.";
		String helpWinSize = "This configures the fragment size m of Aligned Fragment Pairs (AFPs).";
		String helpScoring = "Which scoring function to use.";
		params.add(helpMaxGap);
		//params.add(helpRmsdThr);
		params.add(helpWinSize);
		params.add(helpScoring);
		params.add("The maximum RMSD at which to stop alignment optimization. (default: unlimited=99)");
		params.add("Gap opening penalty during alignment optimization.");
		params.add("Gap extension penalty during alignment optimization.");
		return params;
	}

	public List<String> getUserConfigParameters() {
		List<String> params = new ArrayList<String>();
		params.add("MaxGapSize");
		//params.add("RmsdThr");
		params.add("WinSize");
		params.add(SCORING_STRATEGY);
		params.add("MaxOptRMSD");
		params.add("GapOpen");
		params.add("GapExtension");

		return params;
	}

	public List<String> getUserConfigParameterNames(){
		List<String> params = new ArrayList<String>();
		params.add("max. gap size G (during AFP extension).");
		//params.add("RMSD threshold during trace of the fragment matrix.");
		params.add("fragment size m");
		params.add("Which scoring function to use");
		params.add("RMSD threshold for alignment.");
		params.add("Gap open");
		params.add("Gap extension");
		return params;
	}

	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = new ArrayList<Class>();
		params.add(Integer.class);
		//params.add(Double.class);
		params.add(Integer.class);
		params.add(Integer.class);
		params.add(Double.class);
		params.add(Double.class);
		params.add(Double.class);
		return params;
	}



	/**
	 * @return whether information about AFPs should be printed
	 */
	public boolean isShowAFPRanges()
	{
		return showAFPRanges;
	}
	public void setShowAFPRanges(boolean showAFPRanges)
	{
		this.showAFPRanges = showAFPRanges;
	}





	/** set the maximum RMSD cutoff to be applied during alignment optimization. (default: 99 = unlimited)
	 * 
	 * @param param maxOptRMSD
	 */
	public void setMaxOptRMSD(Double param){
		if ( param == null)
			param = 99d;
		maxOptRMSD = param;
	}

	/** Returns the maximum RMSD cutoff to be applied during alignment optimization (default: 99 = unlimited)
	 * 
	 * @return maxOptRMSD
	 */
	public Double getMaxOptRMSD()
	{
		return maxOptRMSD;
	}



	public Double getGapOpen()
	{
		return gapOpen;
	}



	public void setGapOpen(Double gapOpen)
	{
		this.gapOpen = gapOpen;
	}



	public Double getGapExtension()
	{
		return gapExtension;
	}



	public void setGapExtension(Double gapExtension)
	{
		this.gapExtension = gapExtension;
	}



	public Double getDistanceIncrement()
	{
		return distanceIncrement;
	}



	public void setDistanceIncrement(Double distanceIncrement)
	{
		this.distanceIncrement = distanceIncrement;
	}



	/** Get the Original RMSD threshold from which the alignment optimization is started
	 * 
	 * @return oRMSDThreshold
	 */
	public Double getORmsdThr()
	{
		return oRmsdThr;
	}



	/** Set the Original RMSD threshold from which the alignment optimization is started
	 * 
	 * @param oRmsdThr the threshold
	 */
	public void setORmsdThr(Double oRmsdThr)
	{
		this.oRmsdThr = oRmsdThr;
	}


	/** Get the maximum nr of times the (slow) optimiziation of alignment should iterate. Default: unlimited
	 * 
	 * @param maxNrIterationsForOptimization
	 */
	public int getMaxNrIterationsForOptimization() {
		return maxNrIterationsForOptimization;
	}


	/** Set the maximum nr of times the (slow) optimiziation of alignment should iterate. Default: unlimited
	 * 
	 * @param maxNrIterationsForOptimization
	 */
	public void setMaxNrIterationsForOptimization(int maxNrIterationsForOptimization) {
		this.maxNrIterationsForOptimization = maxNrIterationsForOptimization;
	}


	/** Should sequence conservation be considered as part of the alignment? If yes, this weight factor allows to determine how much.
	 *  By default this is set to 0, meaning no contribution of the sequence alignment score.
	 * 
	 * @return seqWeight the weight factor (default 0)
	 */

	public double getSeqWeight() {
		return seqWeight;
	}


	/** Should sequence conservation be considered as part of the alignment? If yes, this weight factor allows to determine how much.
	 *  By default this is set to 0, meaning no contribution of the sequence alignment score.
	 * 
	 * @param seqWeight the weight factor (default 0)
	 */
	public void setSeqWeight(double seqWeight) {
		this.seqWeight = seqWeight;
	}


	/** Sets the  substitution matrix to be used for influencing the alignment with sequence conservation information.
	 * Default: SDM matrix (Prlic et al 2000)
	 * @return substitutionMatrix 
	 */
	public SubstitutionMatrix<AminoAcidCompound> getSubstitutionMatrix() {
		if ( substitutionMatrix == null){
			String matrixName = "PRLA000101";
			substitutionMatrix = SubstitutionMatrixHelper.getMatrixFromAAINDEX(matrixName);

		}
		return substitutionMatrix;
	}


	/** Sets the  substitution matrix to be used for influencing the alignment with sequence conservation information.
	 * Default: SDM matrix (Prlic et al 2000)
	 * @param substitutionMatrix 
	 */
	public void setSubstitutionMatrix(
			SubstitutionMatrix<AminoAcidCompound> substitutionMatrix) {
		this.substitutionMatrix = substitutionMatrix;
	}




}
