/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * http://www.ncbi.nlm.nih.gov/pubmed/14534198
 * </pre>
 * 
 * Thanks to Yuzhen Ye and A. Godzik for granting permission to freely use and redistribute this code.
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
 * Created by Andreas Prlic - RCSB PDB 
 * 
 */

package org.biojava.bio.structure.align.fatcat.calc;

import java.io.StringWriter;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;


public class FatCatParameters implements ConfigStrucAligParams
{

	public static final int DEFAULT_FRAGLEN = 8;

	int fragLen  ; // the length of the fragments to consider...
	int fragLenSq ;
	Double rmsdCut; // cutoff for AFP detection.
	double disCut; // for AFPs connection, to be tuned, 4.0
	double afpDisCut;
	double afpDisCut0;
	double disSmooth; // for smooth calculation of twist penalty calculation
	int misCut;
	int maxGap;
	int maxGapFrag;
	double disFilter;
	double badRmsd;
	int maxTra; // the maximum number of Twists that are allowed...
	double gapCreate;
	double gapExtend;
	double misScore;
	double torsionPenalty;
	double maxPenalty;
	double resScore;
	double fragScore;
	int sparse;

	public FatCatParameters(){
		reset();
	}


	public void reset(){
		fragLen = DEFAULT_FRAGLEN;
		fragLenSq = fragLen * fragLen;
		rmsdCut = 3.0; //cutoff for AFP detection
		disCut = 5.0; //for AFPs connection, to be tuned, 4.0
		afpDisCut = fragLenSq * disCut * disCut;
		afpDisCut0 = fragLenSq * disCut;
		disSmooth = 4.0; //for smoothly calculation of twist penalty calculation
		misCut = 2 * fragLen; //structural-dismilar ranges allowed between AFPs
		maxGap = 40; //try-1 30
		maxGapFrag = fragLen + maxGap;
		disFilter = 2.0 * rmsdCut; //for single AFP denifition to be tuned! //two CA-dis is 3.6
		badRmsd = 4.0; //very important paramerter for twists detection
		maxTra = 5;
		gapCreate = -5.0;
		gapExtend = -0.5;
		misScore = gapExtend; //comparable to gapExtend
		torsionPenalty = 5 * gapCreate; //to be tuned
		maxPenalty = 1 * gapCreate; //to be tuned
		resScore = 3.0; //on average, the score for each well-matched residue pair
		fragScore = resScore * fragLen; //the score for each well-matched fragment
		sparse = 0;
	}


	public Integer getFragLen()
	{
		return fragLen;
	}


	public void setFragLen(Integer fragLen)
	{
		this.fragLen = fragLen;
	}


	public int getFragLenSq()
	{
		return fragLenSq;
	}


	public void setFragLenSq(int fragLenSq)
	{
		this.fragLenSq = fragLenSq;
	}


	/** The cutoff to be used during AFP detection
	 * 
	 * @return rmsdCut parameter
	 */
	public Double getRmsdCut()
	{
		return rmsdCut;
	}

	/** The cutoff to be used during AFP detection
	 * 
	 * @param rmsdCut
	 */
	public void setRmsdCut(Double rmsdCut)
	{
		this.rmsdCut = rmsdCut;
	}

	/** Get the distance cutoff used during AFP chain connectivity checks
	 * 
	 * @return distance Cutoff
	 */
	public Double getDisCut()
	{
		return disCut;
	}


	public void setDisCut(Double disCut)
	{
		this.disCut = disCut;
	}


	public double getAfpDisCut()
	{
		return afpDisCut;
	}


	public void setAfpDisCut(double afpDisCut)
	{
		this.afpDisCut = afpDisCut;
	}


	public double getAfpDisCut0()
	{
		return afpDisCut0;
	}


	public void setAfpDisCut0(double afpDisCut0)
	{
		this.afpDisCut0 = afpDisCut0;
	}


	public double getDisSmooth()
	{
		return disSmooth;
	}


	public void setDisSmooth(double disSmooth)
	{
		this.disSmooth = disSmooth;
	}


	public int getMisCut()
	{
		return misCut;
	}


	public void setMisCut(int misCut)
	{
		this.misCut = misCut;
	}


	public int getMaxGap()
	{
		return maxGap;
	}


	public void setMaxGap(int maxGap)
	{
		this.maxGap = maxGap;
	}


	public int getMaxGapFrag()
	{
		return maxGapFrag;
	}


	public void setMaxGapFrag(int maxGapFrag)
	{
		this.maxGapFrag = maxGapFrag;
	}


	public double getDisFilter()
	{
		return disFilter;
	}


	public void setDisFilter(double disFilter)
	{
		this.disFilter = disFilter;
	}


	public double getBadRmsd()
	{
		return badRmsd;
	}


	public void setBadRmsd(double badRmsd)
	{
		this.badRmsd = badRmsd;
	}


	/** get the maximum number of Twists that are allowed...
	 * 
	 * @return max nr of allowed twists
	 */
	 public Integer getMaxTra()
	{
		return maxTra;
	}

	/** set the maximum number of Twists that are allowed...
	 * 
	 * @param maxTra
	 */
	 public void setMaxTra(Integer maxTra)
	 {
		 this.maxTra = maxTra;
	 }


	 public double getGapCreate()
	 {
		 return gapCreate;
	 }


	 public void setGapCreate(double gapCreate)
	 {
		 this.gapCreate = gapCreate;
	 }


	 public double getGapExtend()
	 {
		 return gapExtend;
	 }


	 public void setGapExtend(double gapExtend)
	 {
		 this.gapExtend = gapExtend;
	 }


	 public double getMisScore()
	 {
		 return misScore;
	 }


	 public void setMisScore(double misScore)
	 {
		 this.misScore = misScore;
	 }


	 public double getTorsionPenalty()
	 {
		 return torsionPenalty;
	 }


	 public void setTorsionPenalty(double torsionPenalty)
	 {
		 this.torsionPenalty = torsionPenalty;
	 }


	 public double getMaxPenalty()
	 {
		 return maxPenalty;
	 }


	 public void setMaxPenalty(double maxPenalty)
	 {
		 this.maxPenalty = maxPenalty;
	 }


	 public double getResScore()
	 {
		 return resScore;
	 }


	 public void setResScore(double resScore)
	 {
		 this.resScore = resScore;
	 }


	 public double getFragScore()
	 {
		 return fragScore;
	 }


	 public void setFragScore(double fragScore)
	 {
		 this.fragScore = fragScore;
	 }


	 public int getSparse()
	 {
		 return sparse;
	 }


	 public void setSparse(int sparse)
	 {
		 this.sparse = sparse;
	 }


	 public List<String> getUserConfigHelp() {
		 List<String> params = new ArrayList<String>();
		 String fragLen = "The length of the fragments.";
		 String rmsdCutHelp = "The RMSD cutoff to be used during AFP detection.";
		 String disCutHelp = "The distance cutoff used when calculate the connectivity of AFP pairs";
		 String twistHelp ="The number of twists that are allowed to be introduced. If set to 0 alignments are run in RIGID mode.";
		 params.add(fragLen);
		 params.add(rmsdCutHelp);
		 params.add(disCutHelp);
		 params.add(twistHelp);
		 return params;

	 }


	 public List<String> getUserConfigParameterNames() {
		 List<String> params = new ArrayList<String>();
		 params.add("Fragment Length");
		 params.add("RMSD Cutoff");
		 params.add("AFP Distance Cutoff");
		 params.add("Maximum Nr. of twists");
		 return params;
	 }


	 public List<String> getUserConfigParameters() {
		 List<String> params = new ArrayList<String>();
		 params.add("FragLen");
		 params.add("RmsdCut");
		 params.add("DisCut");
		 params.add("MaxTra");
		 return params;
	 }


	 @SuppressWarnings({  "rawtypes" })
	public List<Class> getUserConfigTypes() {
		 
		 List<Class> params = new ArrayList<Class>();
		 params.add(Integer.class);
		 params.add(Double.class);
		 params.add(Double.class);
		 params.add(Integer.class);
		 return params;
	 }

	 
	 public String toString(){
	    StringWriter writer = new StringWriter();
	    writer.append("[");
	    if ( maxTra == 0)
	       writer.append("Mode: rigid, ");
	    else 
	       writer.append("Mode: flexible, ");
	    List<String> params = getUserConfigParameters();
	    
	    for ( String s : params){
	       writer.append(s);
	       writer.append(": ");
	       Object val = getValue(s);
	       writer.append(val.toString());
	       writer.append(", ");
	    }
	    writer.append("]");
	    return writer.toString();
	 }
	 
	 @SuppressWarnings({ "unchecked", "rawtypes" })
	   private Object  getValue(String name){

	      try {
	         String methodName = "get" + name;

	         Class paramC = this.getClass();

	         Method m =paramC.getMethod(methodName,(Class[])null);

	         Object value = m.invoke(this);

	         return value;
	      } catch (Exception e){
	         e.printStackTrace();
	         return null;
	      }


	   }

}
