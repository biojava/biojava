/*
 *                    PDB web development code
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
 *
 * Created on Jun 17, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.jama.Matrix;


/** a bean to contain the core of an alignment
 * 
 * The FatCat aligner class is working on the AFPChain class.
 * 
 * @author Andreas Prlic
 * 
 *
 */

public class AFPChain implements Serializable, Cloneable
{

	/**
	 * 
	 */
	private static final long serialVersionUID = -4474029015606617947L;

	public static final String newline = System.getProperty("line.separator");

	/** the default algorithm used in the RCSB PDB all vs. all database searches
	 * 
	 */
	public static final String DEFAULT_ALGORITHM_NAME = "jFatCat_rigid";

	String algorithmName;

	String version;

	String name1;
	String name2;
	long ioTime;
	long calculationTime;
	long id;
	// results:
	double alignScore;
	double alignScoreUpdate;
	int afpChainTwiNum;
	// end of results
	double tmScore;
	// utility
	int minLen ; // the length of the shorter 2 proteins.


	List<AFP> afpSet;
	int[][] afpIndex;
	int[][] afpAftIndex;
	int[][] afpBefIndex;

	Matrix disTable1;
	Matrix disTable2;

	int[] twi = null ; //the number of twists making the best score ending at each AFP

	int afpChainLen;
	int[] afpChainList;
	double[] afpChainTwiBin;
	double[] afpChainTwiList;
	double chainRmsd;

	int chainLen,misLen,gapLen;
	int     blockNum;       //the final block number
	int     blockNumIni;    //block number before block clustering and split
	int     blockNumClu;    //block number after clustering blocks
	int     blockNumSpt;    //block number after spliting blocks
	double[]  blockRmsd;     //the RMSD of each block
	int[]     block2Afp;     //the index of afp for each block
	int[]     blockSize;     //the number of AFPs involved in a block
	double[]  blockScore;    //the score associated with each block
	int[]     blockGap;      //the gaps in each block
	int[]     blockResSize;  //the number of residues involved in a block
	int[][][]     blockResList;//the list of AFP for each block
	Matrix[] blockRotationMatrix;
	Atom[]   blockShiftVector;

	int     focusResn;      //the size of the set
	int[]     focusRes1;     //the residues from protein 1
	int[]     focusRes2;     //the residues from protein 2
	int     focusAfpn;      //the AFP number
	int[]     focusAfpList;  //the AFP list

	boolean shortAlign;

	String [][][] pdbAln; // only needed temp. during XML serialization, since we don;t have coordinates loaded at that time and can map from PDB res numbers to atom positions.
	int[][][] optAln;
	int[] optLen ;
	double[] optRmsd ;
	int optLength;

	char[] alnsymb;
	char[] alnseq1;
	char[] alnseq2;
	int alnLength;
	int alnbeg1;
	int alnbeg2;

	int totalLenIni;
	int totalLenOpt = 0;

	double totalRmsdIni;
	double totalRmsdOpt;

	int ca1Length;
	int ca2Length;

	// this one is special. it comes from the FatCatParameters...
	// default is flexible alignment...
	int maxTra = 5;

	Double conn;
	Double dvar;

	double probability;
	double identity;
	double similarity;
	double normAlignScore;

	int myResultsEQR;
	int myResultsSimilarity1;
	int myResultsSimilarity2;

	// Mark whether the alignment topology is sequential
	// false if a circular permutation has occured
	boolean sequentialAlignment;

	// background distances
	Matrix distanceMatrix;
	
	String description2;

	public AFPChain(){
		init();
	}

	/**
	 * Copy constructor
	 * @param o AFPChain to duplicate
	 */
	public AFPChain(AFPChain o) {
		this.algorithmName = o.algorithmName;
		this.version = o.version;
		this.name1 = o.name1;
		this.name2 = o.name2;
		this.ioTime = o.ioTime;
		this.calculationTime = o.calculationTime;
		this.id = o.id;
		this.alignScore = o.alignScore;
		this.alignScoreUpdate = o.alignScoreUpdate;
		this.afpChainTwiNum = o.afpChainTwiNum;
		this.minLen = o.minLen;
		this.afpSet = new ArrayList<AFP>(o.afpSet);
		this.afpIndex = o.afpIndex == null? null: o.afpIndex.clone();
		this.afpAftIndex = o.afpAftIndex == null? null: o.afpAftIndex.clone();
		this.afpBefIndex = o.afpBefIndex == null? null: o.afpBefIndex.clone();
		this.disTable1 = o.disTable1 == null? null: (Matrix) o.disTable1.clone();
		this.disTable2 = o.disTable2 == null? null: (Matrix) o.disTable2.clone();
		this.twi = o.twi == null ? null : o.twi.clone();
		this.afpChainLen = o.afpChainLen;
		this.afpChainList = o.afpChainList == null? null: o.afpChainList.clone();
		this.afpChainTwiBin = o.afpChainTwiBin == null? null: o.afpChainTwiBin.clone();
		this.afpChainTwiList = o.afpChainTwiList == null? null: o.afpChainTwiList.clone();
		this.chainRmsd = o.chainRmsd;
		this.chainLen = o.chainLen;
		this.misLen = o.misLen;
		this.gapLen = o.gapLen;
		this.blockNum = o.blockNum;
		this.blockNumIni = o.blockNumIni;
		this.blockNumClu = o.blockNumClu;
		this.blockNumSpt = o.blockNumSpt;
		this.blockRmsd = o.blockRmsd == null ? null : o.blockRmsd.clone();
		this.block2Afp = o.block2Afp == null ? null : o.block2Afp.clone();
		this.blockSize = o.blockSize == null ? null : o.blockSize.clone();
		this.blockScore = o.blockScore == null ? null : o.blockScore.clone();
		this.blockGap = o.blockGap == null ? null : o.blockGap.clone();
		this.blockResSize = o.blockResSize == null ? null : o.blockResSize.clone();
		this.blockResList = o.blockResList == null ? null : o.blockResList.clone();
		this.blockRotationMatrix = o.blockRotationMatrix == null ? null : o.blockRotationMatrix.clone();
		this.blockShiftVector = o.blockShiftVector == null ? null : o.blockShiftVector.clone();
		this.focusResn = o.focusResn;
		this.focusRes1 = o.focusRes1 == null ? null : o.focusRes1.clone();
		this.focusRes2 = o.focusRes2 == null ? null : o.focusRes2.clone();
		this.focusAfpn = o.focusAfpn;
		this.focusAfpList = o.focusAfpList == null ? null : o.focusAfpList.clone();
		this.shortAlign = o.shortAlign;
		this.pdbAln = o.pdbAln == null ? null : o.pdbAln.clone();
		this.optAln = o.optAln == null ? null : o.optAln.clone();
		this.optLen = o.optLen == null ? null : o.optLen.clone();
		this.optRmsd = o.optRmsd == null ? null : o.optRmsd.clone();
		this.optLength = o.optLength;
		this.alnsymb = o.alnsymb == null ? null : o.alnsymb.clone();
		this.alnseq1 = o.alnseq1 == null ? null : o.alnseq1.clone();
		this.alnseq2 = o.alnseq2 == null ? null : o.alnseq2.clone();
		this.alnLength = o.alnLength;
		this.alnbeg1 = o.alnbeg1;
		this.alnbeg2 = o.alnbeg2;
		this.totalLenIni = o.totalLenIni;
		this.totalLenOpt = o.totalLenOpt;
		this.totalRmsdIni = o.totalRmsdIni;
		this.totalRmsdOpt = o.totalRmsdOpt;
		this.ca1Length = o.ca1Length;
		this.ca2Length = o.ca2Length;
		this.maxTra = o.maxTra;
		this.conn = new Double(o.conn);
		this.dvar = new Double(o.dvar);
		this.probability = o.probability;
		this.identity = o.identity;
		this.similarity = o.similarity;
		this.normAlignScore = o.normAlignScore;
		this.myResultsEQR = o.myResultsEQR;
		this.myResultsSimilarity1 = o.myResultsSimilarity1;
		this.myResultsSimilarity2 = o.myResultsSimilarity2;
		this.distanceMatrix = o.distanceMatrix;
	}

	/**
	 * Creates and returns a copy of this object.
	 */
	public Object clone() {
		return new AFPChain(this);
	}

	public long getId()
	{
		return id;
	}
	public void setId(long id)
	{
		this.id = id;
	}

	public String toCE(Atom[] ca1, Atom[]ca2) {

		return AfpChainWriter.toCE(this,ca1,ca2);

	}

	public String toRotMat(){

		return AfpChainWriter.toRotMat(this);
	}

	public String toFatcat(Atom[] ca1, Atom[]ca2){

		return AfpChainWriter.toFatCat(this,ca1,ca2);

	}

	public String toDBSearchResult(){

		return AfpChainWriter.toDBSearchResult(this);
	}

	protected void calcSimilarity() {
		Map<String,Double> idMap = AFPAlignmentDisplay.calcIdSimilarity(alnseq1,alnseq2,alnLength);

		//probability = idMap.get("probability");
		similarity = idMap.get("similarity");
		identity   = idMap.get("identity");
		
	}




	/** 
	 * Get the number of structurally equivalent residues
	 * 
	 * @return nr of EQR
	 */
	public int getNrEQR(){

		if (myResultsEQR < 0){
			if ( optLen == null) {
				myResultsEQR = 0;
				return 0;
			}

			int nrEqr = 0;
			for(int bk = 0; bk < blockNum; bk ++)       {        

				for ( int i=0;i< optLen[bk];i++){
					nrEqr++;
				}
			}
			myResultsEQR = nrEqr;
		}
		return myResultsEQR;
	}

	/** Get the coverage of protein 1 with the alignment
	 * 
	 * @return percentage of coverage, between 0 and 100.
	 */
	public int getCoverage1(){
		if ( myResultsSimilarity1 < 0 ) {
			int distance = ca1Length + ca2Length - 2 * getNrEQR();

			int similarity = (ca1Length + ca2Length - distance ) / 2;

			myResultsSimilarity1 = Math.round(similarity /(float) ca1Length * 100);
		} 
		return myResultsSimilarity1;
	}
	
	/** Get the coverage of protein 2 with the alignment
	 * 
	 * @return percentage of coverage, between 0 and 100.
	*/
	public int getCoverage2(){
		if ( myResultsSimilarity2 < 0 ) {

			
			int distance = ca1Length + ca2Length - 2 * getNrEQR();

			int similarity = (ca1Length + ca2Length - distance ) / 2;
			myResultsSimilarity2 = Math.round(similarity /(float) ca2Length * 100);
		}
		return myResultsSimilarity2;

	}
	
	/** get the coverage of protein 1 with the alignment
	 * 
	 * @return percentage of coverage, between 0 and 100.
	 * @deprecated use getCoverage1() instead
	 */
	@Deprecated 
	public int getSimilarity1(){
		return getCoverage1();
		
	}

	/** get the coverage of protein 2 with the alignment
	 * 
	 * @return percentage of coverage, between 0 and 100.
	 * @deprecated use getCoverage2() instead
	 */
	@Deprecated
	public int getSimilarity2(){
		return getCoverage2();
		
	}

	public String toString(){

		//int lA = ca1Length;
		//int lB = ca2Length;
		//int distance = lA + lB - 2 * getNrEQR();

		StringBuffer str = new StringBuffer("");
		str.append("EQR:");
		str.append(getNrEQR());

		str.append("\tLen1:");
		str.append(this.getCa1Length());
		str.append("\tLen2:");
		str.append(this.getCa2Length());
		str.append(String.format("\tscore: %.2f",this.getAlignScore()));
		str.append("\t");		
		if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName) || algorithmName.equalsIgnoreCase(CeSideChainMain.algorithmName)){
			str.append("Z-score:");
			str.append(String.format("%.2f",this.getProbability()));
		} else {
			str.append("Probability:");
			str.append(String.format("%.2e",this.getProbability()));
		}
		str.append("\tRMSD:");
		str.append(String.format("%.2f",this.getTotalRmsdOpt()));

		str.append("\tSeqID:");
        str.append(String.format("%.0f",getIdentity()*100));
		str.append("%\tSeqSim:");
		str.append(String.format("%.0f",getSimilarity()*100));
		str.append("%\tCov1:");
		str.append(this.getCoverage1());
		str.append("%\tCov2:");
		str.append(this.getCoverage2());
		str.append("%");
		
		if (  tmScore != -1)  {
			str.append("\ttmScore:");
			str.append(String.format("%.2f",tmScore));
		}
		str.append(newline);


		return str.toString();
	}

	public boolean isSignificantResult(){
		if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName) || algorithmName.equalsIgnoreCase(CeSideChainMain.algorithmName)){
			if (probability >= 3.5)
				return true;			
		} else {
			if (probability < 0.01)
				return true;
		}
		return false;
	}




	private void init(){
		shortAlign = false;
		afpChainLen = 0;

		afpChainList      = null;
		afpChainTwiBin  = null;
		afpChainTwiList = null;
		chainRmsd=0;
		chainLen = misLen = gapLen = 0;

		blockResSize = null;
		blockScore = null;
		blockGap = null;
		optAln = null;
		pdbAln = null;
		optLen = null;

		optRmsd = null;

		block2Afp = new int[maxTra+1];
		blockSize = new int[maxTra+1];      
		blockRmsd = new double[maxTra+1];
		blockScore = new double[maxTra+1];
		blockGap = new int[maxTra+1];

		blockResSize = new int[maxTra+1];

		afpSet = new ArrayList<AFP>();
		totalLenIni = totalLenOpt = 0;
		totalRmsdIni = totalRmsdOpt = 0.0;

		afpChainTwiNum = 0;
		alignScore = 0;
		alignScoreUpdate = 0;
		conn = new Double(0);
		dvar = new Double(0);
		calculationTime = 0;
		
		similarity = -1;
		identity   = -1;
		myResultsEQR = -1;
		myResultsSimilarity1 = -1;
		myResultsSimilarity2 = -1;
		algorithmName = DEFAULT_ALGORITHM_NAME ;
		version = "1.0";
		sequentialAlignment = true;
		distanceMatrix = null;
		tmScore = -1;
		description2=null;
	}

	/**
	 * Resets properties which can be calculated on the fly
	 */
	private void invalidate() {
		myResultsEQR = -1;
		myResultsSimilarity1 = -1;
		myResultsSimilarity2 = -1;
		identity = -1;
		similarity = -1;
	}

	/** used temporarily during XML serialization to track the PDB positions of the alignmnet
	 * 
	 * @return String array
	 */
	public String[][][] getPdbAln() {
		return pdbAln;
	}


	public void setPdbAln(String[][][] pdbAln) {
		this.pdbAln = pdbAln;
	}


	public Double getConn()
	{
		return conn;
	}


	public void setConn(Double conn)
	{
		this.conn = conn;
	}


	public Double getDVar()
	{
		return dvar;
	}


	public void setDVar(Double dvar)
	{
		this.dvar = dvar;
	}


	/** get the maximum nr of Twists that are allowed...
	 * 
	 * @return maxTra, the max nr of twists 
	 */
	public int getMaxTra()
	{
		return maxTra;
	}

	/**
	 * Set the maximum number of Twists that are allowed...
	 * @param maxTra
	 */
	public void setMaxTra(int maxTra)
	{
		this.maxTra = maxTra;
	}


	public double getAlignScore()
	{
		return alignScore;
	}

	public void setAlignScore(double alignScore)
	{
		this.alignScore = alignScore;
	}

	public double getAlignScoreUpdate()
	{
		return alignScoreUpdate;
	}

	public void setAlignScoreUpdate(double alignScoreUpdate)
	{
		this.alignScoreUpdate = alignScoreUpdate;
	}

	public int getAfpChainTwiNum()
	{
		return afpChainTwiNum;
	}

	public void setAfpChainTwiNum(int afpChainTwiNum)
	{
		this.afpChainTwiNum = afpChainTwiNum;
	}

	public int getMinLen()
	{
		return minLen;
	}

	public void setMinLen(int minLen)
	{
		this.minLen = minLen;
	}



	/**
	 * Get the set of AFPs for this alignment.
	 * An AFP is a local ungapped alignment between the two peptides.
	 * 
	 * AFPs are set before the final optimization step. To get the final
	 * alignment, look at the aligned pairs from {@link #getOptAln()}.
	 *  
	 * @return The optimal set of AFPs
	 * @see #getOptAln()
	 */
	public List<AFP> getAfpSet()
	{
		return afpSet;
	}


    /**
     * Set the set of AFPs for this alignment.
     * An AFP is a local ungapped alignment between the two peptides.
     * 
     * AFPs are set before the final optimization step. To get the final
     * alignment, look at the aligned pairs from {@link #getOptAln()}.
     */  
	public void setAfpSet(List<AFP> afpSet)
	{
		this.afpSet = afpSet;
	}

	public int[][] getAfpIndex()
	{
		return afpIndex;
	}

	public void setAfpIndex(int[][] afpIndex)
	{
		this.afpIndex = afpIndex;
	}


	public int[][] getAfpAftIndex()
	{
		return afpAftIndex;
	}

	public void setAfpAftIndex(int[][] afpAftIndex)
	{
		this.afpAftIndex = afpAftIndex;
	}


	public int[][] getAfpBefIndex()
	{
		return afpBefIndex;
	}

	public void setAfpBefIndex(int[][] afpBefIndex)
	{
		this.afpBefIndex = afpBefIndex;
	}


	public Matrix getDisTable1()
	{
		return disTable1;
	}

	public void setDisTable1(Matrix disTable1)
	{
		this.disTable1 = disTable1;
	}


	public Matrix getDisTable2()
	{
		return disTable2;
	}

	public void setDisTable2(Matrix disTable2)
	{
		this.disTable2 = disTable2;
	}


	public int[] getTwi()
	{
		return twi;
	}

	public void setTwi(int[] twi)
	{
		this.twi = twi;
	}

	public int getAfpChainLen()
	{
		return afpChainLen;
	}

	public void setAfpChainLen(int afpChainLen)
	{
		this.afpChainLen = afpChainLen;
	}

	public int[] getAfpChainList()
	{
		return afpChainList;
	}

	public void setAfpChainList(int[] afpChainList)
	{
		this.afpChainList = afpChainList;
	}

	public double[] getAfpChainTwiBin()
	{
		return afpChainTwiBin;
	}

	public void setAfpChainTwiBin(double[] afpChainTwiBin)
	{
		this.afpChainTwiBin = afpChainTwiBin;
	}

	public double[] getAfpChainTwiList()
	{
		return afpChainTwiList;
	}

	public void setAfpChainTwiList(double[] afpChainTwiList)
	{
		this.afpChainTwiList = afpChainTwiList;
	}

	public double getChainRmsd()
	{
		return chainRmsd;
	}

	/** The RMSD of the chain of AFPs. Set during AFPCHainer.traceBack();
	 * 
	 * @param chainRmsd
	 */
	public void setChainRmsd(double chainRmsd)
	{
		this.chainRmsd = chainRmsd;
	}

	public int getChainLen()
	{
		return chainLen;
	}

	public void setChainLen(int chainLen)
	{
		this.chainLen = chainLen;
	}

	public int getMisLen()
	{
		return misLen;
	}

	public void setMisLen(int misLen)
	{
		this.misLen = misLen;
	}

	public int getGapLen()
	{
		return gapLen;
	}

	public void setGapLen(int gapLen)
	{
		this.gapLen = gapLen;
	}

	/** The number of blocks in the alignment
	 * 
	 * @return the nr of blocks in alignment
	 */
	public int getBlockNum()
	{
		return blockNum;
	}

	public void setBlockNum(int blockNum)
	{         
		this.blockNum = blockNum;
	}

	public int getBlockNumIni()
	{
		return blockNumIni;
	}

	public void setBlockNumIni(int blockNumIni)
	{
		this.blockNumIni = blockNumIni;
	}

	public int getBlockNumClu()
	{
		return blockNumClu;
	}

	public void setBlockNumClu(int blockNumClu)
	{
		this.blockNumClu = blockNumClu;
	}

	public int getBlockNumSpt()
	{
		return blockNumSpt;
	}

	public void setBlockNumSpt(int blockNumSpt)
	{
		this.blockNumSpt = blockNumSpt;
	}

	public double[] getBlockRmsd()
	{
		return blockRmsd;
	}

	public void setBlockRmsd(double[] blockRmsd)
	{
		this.blockRmsd = blockRmsd;
	}

	public int[] getBlock2Afp()
	{
		return block2Afp;
	}

	public void setBlock2Afp(int[] block2Afp)
	{
		this.block2Afp = block2Afp;
	}

	public int[] getBlockSize()
	{
		return blockSize;
	}

	public void setBlockSize(int[] blockSize)
	{
		this.blockSize = blockSize;
	}

	public double[] getBlockScore()
	{
		return blockScore;
	}

	public void setBlockScore(double[] blockScore)
	{
		this.blockScore = blockScore;
	}

	public int[] getBlockGap()
	{
		return blockGap;
	}

	public void setBlockGap(int[] blockGap)
	{
		this.blockGap = blockGap;
	}

	public int[] getBlockResSize()
	{
		return blockResSize;
	}

	public void setBlockResSize(int[] blockResSize)
	{
		this.blockResSize = blockResSize;
	}


	/** tracks the residues of the initial blocks (before optimization)
	 * 
	 * 
	 * @return list of block residues
	 */
	public int[][][] getBlockResList()
	{
		return blockResList;
	}

	public void setBlockResList(int[][][] blockResList)
	{
		this.blockResList = blockResList;
	}

	public int getFocusResn()
	{
		return focusResn;
	}

	public void setFocusResn(int focusResn)
	{
		this.focusResn = focusResn;
	}


	public int[] getFocusRes1()
	{
		return focusRes1;
	}

	public void setFocusRes1(int[] focusRes1)
	{
		this.focusRes1 = focusRes1;
	}


	public int[] getFocusRes2()
	{
		return focusRes2;
	}

	public void setFocusRes2(int[] focusRes2)
	{
		this.focusRes2 = focusRes2;
	}

	public int getFocusAfpn()
	{
		return focusAfpn;
	}

	public void setFocusAfpn(int focusAfpn)
	{
		this.focusAfpn = focusAfpn;
	}

	public int[] getFocusAfpList()
	{
		return focusAfpList;
	}

	public void setFocusAfpList(int[] focusAfpList)
	{
		this.focusAfpList = focusAfpList;
	}

	public boolean isShortAlign()
	{
		return shortAlign;
	}

	public void setShortAlign(boolean shortAlign)
	{
		this.shortAlign = shortAlign;
	}

	/** Tracks the Atom positions in the optimal alignment. Note: only considers the equivalent positions, gaps are ignored...
	 * first dimension is the block nr
	 * second dimension is 0 or 1 (the alignment chain index)
	 * third is the position
	 * @return int array
	 */
	public int[][][] getOptAln()
	{
		return optAln;
	}

	public void setOptAln(int[][][] optAln)
	{
		invalidate();
		this.optAln = optAln;
	}

	/**
	 * The length of each block
	 * @return lengths
	 */
	public int[] getOptLen()
	{
		return optLen;
	}

	public void setOptLen(int[] optLen)
	{
		this.optLen = optLen;
	}

	public double[] getOptRmsd()
	{
		return optRmsd;
	}

	public void setOptRmsd(double[] optRmsd)
	{
		this.optRmsd = optRmsd;
	}

	public int getOptLength()
	{
		return optLength;
	}

	/** The length of the optimal alignment. Set by AFPOptimizer.optimizeAln().
	 * This should be the sum of the elements in optLen
	 * @param optLength 
	 */
	public void setOptLength(int optLength)
	{
		this.optLength = optLength;
	}


	public char[] getAlnsymb()
	{
		return alnsymb;
	}

	public void setAlnsymb(char[] alnsymb)
	{
		this.alnsymb = alnsymb;
	}


	public char[] getAlnseq1()
	{
		return alnseq1;
	}

	public void setAlnseq1(char[] alnseq1)
	{
		this.alnseq1 = alnseq1;
	}


	public char[] getAlnseq2()
	{
		return alnseq2;
	}

	public void setAlnseq2(char[] alnseq2)
	{
		this.alnseq2 = alnseq2;
	}


	/**
	 * @return The total length of the alignment, including gaps
	 * @see #getOptLength(), the number of aligned residues in the final alignment.
	 */
	public int getAlnLength()
	{
		return alnLength;
	}

	public void setAlnLength(int alnLength)
	{
		this.alnLength = alnLength;
	}

	/**
	 * @return The index of the first aligned residue in protein 1
	 */
	public int getAlnbeg1()
	{
		return alnbeg1;
	}

	public void setAlnbeg1(int alnbeg1)
	{
		this.alnbeg1 = alnbeg1;
	}
	/**
	 * @return The index of the first aligned residue in protein 2
	 */
	public int getAlnbeg2()
	{
		return alnbeg2;
	}

	public void setAlnbeg2(int alnbeg2)
	{
		this.alnbeg2 = alnbeg2;
	}

	public int getTotalLenIni()
	{
		return totalLenIni;
	}

	public void setTotalLenIni(int totalLenIni)
	{
		this.totalLenIni = totalLenIni;
	}

	public int getTotalLenOpt()
	{
		return totalLenOpt;
	}

	public void setTotalLenOpt(int totalLenOpt)
	{
		this.totalLenOpt = totalLenOpt;
	}

	/** this is the init-RMSD, not the final RMSD after refinement.
	 * 
	 * @return totalRmsdIni
	 */
	public double getTotalRmsdIni()
	{
		return totalRmsdIni;
	}

	/** this is the init-RMSD, not the final RMSD after refinement.
	 * 
	 * @param totalRmsdIni
	 */
	public void setTotalRmsdIni(double totalRmsdIni)
	{
		this.totalRmsdIni = totalRmsdIni;
	}


	/** The RMSD of the final alignment. Use this to print overal alignment RMSD.
	 * 
	 * @return total RMSD of the optimal alignment.
	 */
	public double getTotalRmsdOpt()
	{
		return totalRmsdOpt;
	}

	/** The RMSD of the final alignment. Use this to print overal alignment RMSD.
	 * 
	 * @param totalRmsdOpt : total RMSD of the optimal alignment
	 */
	public void setTotalRmsdOpt(double totalRmsdOpt)
	{
		this.totalRmsdOpt = totalRmsdOpt;
	}


	public String getName1()
	{
		return name1;
	}


	public void setName1(String name1)
	{
		this.name1 = name1;
	}



	public String getName2()
	{
		return name2;
	}

	public void setName2(String name2)
	{
		this.name2 = name2;
	}


	public long getCalculationTime()
	{
		return calculationTime;
	}

	public void setCalculationTime(long calculationTime)
	{
		this.calculationTime = calculationTime;
	}

	public int getCa1Length()
	{
		return ca1Length;
	}

	public void setCa1Length(int ca1Length)
	{
		this.ca1Length = ca1Length;
	}

	public int getCa2Length()
	{
		return ca2Length;
	}

	public void setCa2Length(int ca2Length)
	{
		this.ca2Length = ca2Length;
	}

	public long getIoTime()
	{
		return ioTime;
	}

	public void setIoTime(long ioTime)
	{
		this.ioTime = ioTime;
	}

	/** The probability (FATCAT) or Z-score (CE) of the alignment.
	 * 
	 * @return either the probability (FATCAT) or the Z-score (CE) of the alignment.
	 */
	public double getProbability()
	{
		return probability;
	}

	public void setProbability(double probability)
	{
		this.probability = probability;
	}

/** The percent of residues that are sequence-identical in the alignment.
 * 
 * @return a value between 0 and 1
 */
	public double getIdentity() {
		if ( identity <= 0) {
			System.out.println("recaclulating ID and SIM (" + identity +")");
			calcSimilarity();
		}
		return identity;
	}

	public void setIdentity(double identity) {
		
		this.identity = identity;
	}


	/** Returns the similarity score for the alignment. This gives the percent of 
	 * sequence similar residues in the alignment.
	 * 
	 * @return a double between 0 and 1
	 */
	public double getSimilarity() {
		if ( similarity < 0)
			calcSimilarity();
		return similarity;
	}

	public void setSimilarity(double similarity) {
		this.similarity = similarity;
	}


	public double getNormAlignScore()
	{
		return normAlignScore;
	}

	public void setNormAlignScore(double normAlignScore)
	{
		this.normAlignScore = normAlignScore;
	}

	public Matrix[] getBlockRotationMatrix()
	{
		return blockRotationMatrix;
	}

	public void setBlockRotationMatrix(Matrix[] blockRotationMatrix)
	{
		this.blockRotationMatrix = blockRotationMatrix;
	}

	public Atom[] getBlockShiftVector()
	{
		return blockShiftVector;
	}

	public void setBlockShiftVector(Atom[] blockShiftVector)
	{
		this.blockShiftVector = blockShiftVector;
	}

	public String getAlgorithmName() {
		return algorithmName;
	}

	public void setAlgorithmName(String algorithmName) {
		this.algorithmName = algorithmName;
	}

	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		this.version = version;
	}

	/**
	 * Get whether this alignment has the normal topology, ie the residues
	 * aligned in each block increase sequentially over the original protein.
	 * 
	 * This will be false if a circular permutation was detected.
	 * @return true if the alignment is sequential
	 */
	public boolean isSequentialAlignment() {
		return sequentialAlignment;
	}
	/**
	 * Set whether this alignment has the normal topology, ie the residues
	 * aligned in each block increase sequentially over the original protein.
	 * 
	 * This will be false if a circular permutation was detected.
	 */
	public void setSequentialAlignment(boolean sequentialAlignment) {
		this.sequentialAlignment = sequentialAlignment;
	}

	/**
	 * A matrix with <i>ca1length</i> rows and <i>ca2length</i> columns.
	 * For CE this is the distance matrix, but the exact interpretation is left
	 * up to the alignment algorithm.
	 * 
	 * <p>Note:
	 * A {@link org.biojava.bio.structure.gui.JMatrixPanel}, which is used in
	 * the structure-gui package to display distance matrices, will display the
	 * transpose of this matrix. Be sure to take that into account when debugging
	 * visually.
	 * 
	 * @return A matrix with dimensions ca1length x ca2length, or null
	 */
	public Matrix getDistanceMatrix() {
		return distanceMatrix;
	}

	/**
	 * A matrix with <i>ca1length</i> rows and <i>ca2length</i> columns.
	 * For CE this is the distance matrix, but the exact interpretation is left
	 * up to the alignment algorithm.
	 * @param distanceMatrix A matrix with dimensions ca1length x ca2length
	 */
	public void setDistanceMatrix(Matrix distanceMatrix) {
		this.distanceMatrix = distanceMatrix;
		
		//System.out.println("Setting distMatrix "+(distanceMatrix==null?"null":"not null"));
	}

	
	public void setTMScore(double tmScore){
	   this.tmScore = tmScore;
	}
	
	/** Returns the tmScore of the alignment. If the score has not been calcualted yet,
	 * returns -1. To calculate it call AFPChainScorer.getTMScore(afpChain, ca1, ca2);
	 * 
	 * @return -1, if not calculated, or the TM-score, a score between 0 and 1
	 */
   public double getTMScore()
   {
     
      return tmScore;
   }


   
   /** Get a textual description for the protein 2 of the alignment.
    * 
    * @return
    */
	public String getDescription2() {
		return description2;
	}
	
	
	/** Set the textual description for protein 2.
	 * 
	 * @param desc
	 */
	public void setDescription2(String desc){
		this.description2 = desc;
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + blockNum;
		result = prime * result + ca1Length;
		result = prime * result + ca2Length;
		result = prime * result + Arrays.hashCode(optAln);
		result = prime * result + Arrays.hashCode(optLen);
		result = prime * result + optLength;
		return result;
	}

	/**
	 * A week equality metric.
	 * 
	 * Checks if the optAlign is the same, and if the objects being compared
	 * seem to be the same (same names, lengths). Does not check properties
	 * of the alignment such as scores or superposition matrices.
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AFPChain other = (AFPChain) obj;
		if (blockNum != other.blockNum)
			return false;
		if (ca1Length != other.ca1Length)
			return false;
		if (ca2Length != other.ca2Length)
			return false;
		if (!Arrays.deepEquals(optAln, other.optAln))
			return false;
		if (!Arrays.equals(optLen, other.optLen))
			return false;
		if (optLength != other.optLength)
			return false;
		return true;
	}


}
