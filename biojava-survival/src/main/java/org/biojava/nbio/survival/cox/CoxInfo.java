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
 */
package org.biojava.nbio.survival.cox;

import org.biojava.nbio.survival.cox.stats.ChiSq;
import org.biojava.nbio.survival.kaplanmeier.figure.ExpressionFigure;
import org.biojava.nbio.survival.kaplanmeier.figure.KaplanMeierFigure;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;

//import org.biojava.nbio.survival.cox.comparators.SurvivalInfoComparator;

/**
 * Holds the results of a cox analysis where calling dump(), toString() will give an output similar to R
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxInfo {

	private WaldTestInfo waldTestInfo = null;
	String message = "";
	Integer maxIterations = null;
	Double eps = null;
	Double toler = null;
	CoxMethod method;
	private double[][] imat = null; //:the variance matrix at beta=final
	private double[][] naive_imat = null; //the original variance matrix used in residuals
	double[] u = new double[0]; //:score vector
	int iterations = 0; //:actual number of iterations used
	int flag = 0; // success flag  1000  did not converge 1 to nvar: rank of the solution
	double logTest = 0;
	double logTestpval = 0;
	double loglikInit = 0; //loglik at beta=initial values, at beta=final
	double loglikFinal = 0;
	Double scoreLogrankTest;
	private Double rscore = null; //robust score
	private Double rscoreLogrankTestpvalue = null;
	private double degreeFreedom;
	private Double scoreLogrankTestpvalue;
	int numSamples = 0;
	int numEvents = 0;
	private LinkedHashMap<String, String> metaDataFilter = null;
	private LinkedHashMap<String, CoxCoefficient> coefficientsList = new LinkedHashMap<String, CoxCoefficient>();
	LinkedHashMap<Double, Double> baselineSurvivorFunction = new LinkedHashMap<Double, Double>();
	ArrayList<SurvivalInfo> survivalInfoList = new ArrayList<SurvivalInfo>();
	/**
	 *
	 */
	public KaplanMeierFigure kmf = null;
	/**
	 *
	 */
	public ExpressionFigure ef = null;

	/**
	 *
	 * @return
	 */
	public ArrayList<SurvivalInfo> getSurvivalInfoList() {
		return survivalInfoList;
	}

	/**
	 *
	 * @param var
	 * @throws Exception
	 */
	public void setVariance(double[][] var) throws Exception {

		//   if (Math.abs(var[0][1] - var[1][0]) > .0000000000001) { //in the CoxCC correction looks like a precision error keeps these from being equal 10-19
		//       throw new Exception("Expecting diagonal to be equal");
		//   }
		imat = new double[var.length][var[0].length];
		for (int i = 0; i < var.length; i++) {
			for (int j = 0; j < var[0].length; j++) {
				imat[i][j] = var[i][j];
			}
		}
		calcSummaryValues();
	}

	/**
	 *
	 * @return
	 */
	public double[][] getVariance() {
		double[][] var = new double[imat.length][imat[0].length];
		for (int i = 0; i < var.length; i++) {
			for (int j = 0; j < var[0].length; j++) {
				var[i][j] = imat[i][j];
			}
		}

		return var;
	}

	/**
	 *
	 * @param var
	 * @throws Exception
	 */
	public void setNaiveVariance(double[][] var) throws Exception {
//        if (var[0][1] != var[1][0]) {
//            throw new Exception("Expecting diagonal to be equal");
//        }


		naive_imat = new double[var.length][var[0].length];
		for (int i = 0; i < var.length; i++) {
			for (int j = 0; j < var[0].length; j++) {
				naive_imat[i][j] = var[i][j];
			}
		}

		calcSummaryValues();
	}

	/**
	 *
	 * @return
	 */
	public double[][] getNaiveVariance() {
		double[][] var = new double[imat.length][imat[0].length];
		for (int i = 0; i < var.length; i++) {
			for (int j = 0; j < var[0].length; j++) {
				var[i][j] = naive_imat[i][j];
			}
		}

		return var;

	}

	/**
	 *
	 * @param data
	 */
	public void setSurvivalInfoList(ArrayList<SurvivalInfo> data) {
		survivalInfoList = data;
		numSamples = data.size();

		for (SurvivalInfo si : data) {
			if (si.getStatus() == 1) {
				numEvents++;
			}
		}
	}

	/**
	 *
	 * @return
	 */
	public double[] getWeighted() {
		double[] weighted = new double[survivalInfoList.size()];
		int p = 0;
		for (SurvivalInfo si : this.survivalInfoList) {
			weighted[p] = si.getWeight();
			p++;
		}
		return weighted;
	}

	/**
	 *
	 * @return
	 */
	public double[][] getVariableResiduals() {
		ArrayList<String> variables = new ArrayList<String>(coefficientsList.keySet());
		double[][] rr = new double[survivalInfoList.size()][variables.size()];
		int p = 0;
		for (SurvivalInfo si : this.survivalInfoList) {
			int i = 0;
			for (String v : variables) {
				rr[p][i] = si.getResidualVariable(v);
				i++;
			}
			p++;
		}

		return rr;
	}

	/**
	 *
	 * @param rr
	 */
	public void setVariableResiduals(double[][] rr) {
		ArrayList<String> variables = new ArrayList<String>(coefficientsList.keySet());

		int p = 0;
		for (SurvivalInfo si : this.survivalInfoList) {
			int i = 0;
			for (String v : variables) {
				si.setResidualVariable(v, rr[p][i]);
				i++;
			}
			p++;
		}

	}

	/**
	 *
	 * @return
	 */
	public int getNumberCoefficients() {
		return coefficientsList.size();
	}

	/**
	 *
	 * @param name
	 * @return
	 */
	public CoxCoefficient getCoefficient(String name) {
		return coefficientsList.get(name);
	}

	/**
	 *
	 * @param name
	 * @param coefficient
	 */
	public void setCoefficient(String name, CoxCoefficient coefficient) {
		coefficientsList.put(name, coefficient);
	}

	/**
	 *
	 * @param header
	 * @param beginLine
	 * @param beginCell
	 * @param endCell
	 * @param endLine
	 * @return
	 */
	public String getCoefficientText(boolean header, String beginLine, String beginCell, String endCell, String endLine) {
		String o = "";
		if (header) {
			String robust = "";
			if (naive_imat != null) {
				robust = beginCell + "robust se" + endCell;
			}
			o = o + beginLine + beginCell + fmtpl("", 9) + endCell + beginCell + fmtpl("coef", 9) + endCell + beginCell + fmtpl("se(coef)", 9) + endCell + robust + beginCell + fmtpl("z", 9) + endCell + beginCell + fmtpl("p-value", 9) + endCell + beginCell + fmtpl("HR", 9) + endCell + beginCell + fmtpl("lower .95", 9) + endCell + beginCell + fmtpl("upper .95", 9) + endCell + endLine;
		}//Coefficients,Coe,StdErr,HR,p-value,HR Lo 95%,HR Hi 95%

		for (CoxCoefficient coe : coefficientsList.values()) {
			String robust = "";
			String stderror = "";
			if (naive_imat != null) {
				stderror = beginCell + fmt(coe.getRobustStdError(), 5, 9) + endCell;
				robust = beginCell + fmt(coe.getStdError(), 5, 9) + endCell;
			} else {
				stderror = beginCell + fmt(coe.getStdError(), 5, 9) + endCell;
			}
			o = o + beginLine + beginCell + fmtpr(coe.getName(), 9) + endCell + beginCell + fmt(coe.getCoeff(), 5, 9) + stderror + robust + endCell + beginCell + fmt(coe.getZ(), 5, 9) + endCell + beginCell + fmt(coe.getPvalue(), 6, 9) + endCell + beginCell + fmt(coe.getHazardRatio(), 3, 9) + endCell + beginCell + fmt(coe.getHazardRatioLoCI(), 3, 9) + endCell + beginCell + fmt(coe.getHazardRatioHiCI(), 3, 9) + endCell + endLine;
		}
		return o;
	}

	/**
	 *
	 * @param d
	 * @param precision
	 * @param pad
	 * @return
	 */
	public static String fmt(Double d, int precision, int pad) {
		if(d == null)
			return "";
		if(Double.isNaN(d))
			return "";
		String value = "";
		DecimalFormat dfe = new DecimalFormat("0.00E0");
		String dpad = "0.";
		double p = 1.0;
		for (int i = 0; i < (precision); i++) {
			dpad = dpad + "0";
			p = p / 10.0;
		}
		DecimalFormat df = new DecimalFormat(dpad);
		if (Math.abs(d) >= p) {
			value = df.format(d);
		} else {
			value = dfe.format(d);
		}
		int length = value.length();
		int extra = pad - length;
		if (extra > 0) {
			for (int i = 0; i < extra; i++) {
				value = " " + value;
			}
		}
		return value;
	}

	/**
	 *
	 */
	private void calcSummaryValues() {

		//beta

		ArrayList<String> variables = new ArrayList<String>(coefficientsList.keySet());
		for (int i = 0; i < variables.size(); i++) {
			String variable = variables.get(i);
			CoxCoefficient coe = coefficientsList.get(variable);
			coe.setStdError(Math.sqrt(imat[i][i])); //values can be updated to reflect new error
			if (naive_imat != null) {
				coe.setRobustStdError(Math.sqrt(naive_imat[i][i]));
			}
			coe.setZ(coe.getCoeff() / coe.getStdError());
			coe.setPvalue(ChiSq.norm(Math.abs(coe.getCoeff() / coe.getStdError())));
			//z <- qnorm((1 + conf.int)/2, 0, 1)
			double z = 1.959964;
			coe.setHazardRatioLoCI(Math.exp(coe.getCoeff() - z * coe.getStdError()));
			coe.setHazardRatioHiCI(Math.exp(coe.getCoeff() + z * coe.getStdError()));
		}

		logTest = -2 * (loglikInit - loglikFinal);
		logTestpval = ChiSq.chiSq(logTest, (int) degreeFreedom);

		scoreLogrankTestpvalue = ChiSq.chiSq(scoreLogrankTest, (int) degreeFreedom);
		if (rscore != null) {
			rscoreLogrankTestpvalue = ChiSq.chiSq(rscore, (int) degreeFreedom);
		}
	}

	/**
	 *
	 */
	public void dump() {

		//need an ordered list for comparing to R dumps

//        ArrayList<SurvivalInfo> orderedSurvivalInfoList = new ArrayList<SurvivalInfo>(survivalInfoList);
//        SurvivalInfoComparator sicSort = new SurvivalInfoComparator();
		//       Collections.sort(orderedSurvivalInfoList,sicSort);


		System.out.println();
		System.out.println("$coef");
		for (CoxCoefficient coe : coefficientsList.values()) {
			System.out.print(coe.getCoeff() + " ");
		}
		System.out.println();
		System.out.println("$means");

		for (CoxCoefficient coe : coefficientsList.values()) {
			System.out.print(coe.getMean() + " ");
		}
		System.out.println();
		System.out.println("$u");

		for (double d : u) {
			System.out.print(d + " ");
		}

		System.out.println();
		System.out.println("$imat");
		for (int i = 0; i < imat.length; i++) {
			for (int j = 0; j < imat[0].length; j++) {
				System.out.print(imat[i][j] + " ");
			}
			System.out.println();
		}

		if (this.naive_imat != null) {
			System.out.println("$naive_imat");
			for (int i = 0; i < naive_imat.length; i++) {
				for (int j = 0; j < naive_imat[0].length; j++) {
					System.out.print(naive_imat[i][j] + " ");
				}
				System.out.println();
			}
		}

		System.out.println();
		System.out.println("$loglik");

		System.out.println(loglikInit + " " + loglikFinal);

		System.out.println();
		System.out.println("$sctest");

		System.out.println(this.scoreLogrankTest);

		System.out.println("$iter");
		System.out.println(this.iterations);

		System.out.println("$flag");
		System.out.println(flag);

		System.out.println();
//        if (false) {
//            System.out.println("ID      LP       Score      Residuals");
//            for (SurvivalInfo si : orderedSurvivalInfoList) {
//                System.out.println(si.getOrder() + " " + si.getLinearPredictor() + " " + si.getScore() + " " + si.getResidual());
//
//            }
//            System.out.println();
//            ArrayList<String> variables = new ArrayList<String>(coefficientsList.keySet());
//            System.out.print("Sample");
//            for (String v : variables) {
//                System.out.print("    " + v);
//            }
//            System.out.println("rr");
//            for (SurvivalInfo si : orderedSurvivalInfoList) {
//                System.out.print(si.getOrder());
//                for (String v : variables) {
//                    System.out.print("   " + si.getResidualVariable(v));
//                }
//                System.out.println();
//            }
//        }

	}

	/**
	 *
	 * @param d
	 * @param pad
	 * @return
	 */
	public String fmtpr(String d, int pad) {
		int length = d.length();
		int extra = pad - length;
		if (extra < 0) {
			extra = 0;
		}
		String v = d;
		for (int i = 0; i < extra; i++) {
			v = v + " ";
		}
		return v;
	}

	/**
	 * Pad left a string with spaces
	 *
	 * @param d
	 * @param pad
	 * @return
	 */
	public String fmtpl(String d, int pad) {
		int length = d.length();
		int extra = pad - length;
		if (extra < 0) {
			extra = 0;
		}
		String v = d;
		for (int i = 0; i < extra; i++) {
			v = " " + v;
		}
		return v;
	}

	@Override
	public String toString() {
		return toString("", " ", "\r\n");
	}

	/**
	 *
	 * @param beginLine
	 * @param del
	 * @param endLine
	 * @return
	 */
	public String toString(String beginLine, String del, String endLine) {



		String o = beginLine + fmtpl("", 9) + fmtpl("Avg", 9) + fmtpl("SD", 9) + endLine;
		for (CoxCoefficient coe : coefficientsList.values()) {
			o = o + beginLine + fmtpr(coe.getName(), 9) + fmt(coe.getMean(), 4, 9) + fmt(coe.getStandardDeviation(), 4, 9) + endLine;

		}

		o = o + beginLine + endLine;



		o = o + beginLine + "n= " + this.numSamples + ", number of events=" + this.numEvents + endLine;
		o = o + getCoefficientText(true, beginLine, del, "", endLine);

		o = o + beginLine + endLine;

		if (baselineSurvivorFunction.size() > 0) {
			o = o + beginLine + "Baseline Survivor Function (at predictor means)" + endLine;
			for (Double time : baselineSurvivorFunction.keySet()) {
				Double mean = baselineSurvivorFunction.get(time);
				o = o + beginLine + fmt(time, 4, 10) + fmt(mean, 4, 10) + endLine;
			}
		}

		o = o + beginLine + endLine;
		o = o + beginLine + "Overall Model Fit" + endLine;
		o = o + beginLine + "Iterations=" + iterations + endLine;


		o = o + beginLine + "Likelihood ratio test = " + fmt(this.logTest, 2, 0) + " df=" + this.degreeFreedom + " p-value=" + fmt(this.logTestpval, 7, 0) + endLine;

		o = o + beginLine + "Wald test             = " + fmt(waldTestInfo.getTest(), 2, 0) + " df=" + waldTestInfo.getDf() + " p-value=" + fmt(waldTestInfo.getPvalue(), 7, 0) + endLine;
		o = o + beginLine + "Score (logrank) test  = " + fmt(scoreLogrankTest, 2, 0) + " df=" + ((int) (this.degreeFreedom)) + " p-value=" + fmt(this.scoreLogrankTestpvalue, 7, 0);

		if (this.rscore != null) {
			o = o + ",  Robust = " + fmt(rscore, 2, 0) + " p-value=" + fmt(rscoreLogrankTestpvalue, 7, 0);

		}

		o = o + endLine;

		//       o = o + "Rank of solution flag=" + flag + "\r\n";
		//       o = o + "Log lik Initial=" + loglikInit + "\r\n";
		//       o = o + "Log lik Final=" + loglikFinal + "\r\n";
		o = o + beginLine + "Method=" + method.name() + endLine;




		return o;
	}

	/**
	 * @return the scoreLogrankTest
	 */
	public double getChiSquare() {
		return scoreLogrankTest;
	}

	/**
	 * @return the degreeFreedom
	 */
	public double getDegreeFreedom() {
		return degreeFreedom;
	}

	/**
	 * @return the scoreLogrankTestpvalue
	 */
	public double getOverallModelFitPvalue() {
		return scoreLogrankTestpvalue;
	}

	/**
	 * @return the rscore
	 */
	public Double getRscore() {
		return rscore;
	}

	/**
	 * @param rscore the rscore to set
	 */
	public void setRscore(Double rscore) {
		this.rscore = rscore;
		if (rscore != null) {
			rscoreLogrankTestpvalue = ChiSq.chiSq(rscore, (int) degreeFreedom);
		}
	}

	/**
	 * @return the rscoreLogrankTestpvalue
	 */
	public Double getRscoreLogrankTestpvalue() {
		return rscoreLogrankTestpvalue;
	}

	/**
	 * @param rscoreLogrankTestpvalue the rscoreLogrankTestpvalue to set
	 */
	public void setRscoreLogrankTestpvalue(Double rscoreLogrankTestpvalue) {
		this.rscoreLogrankTestpvalue = rscoreLogrankTestpvalue;
	}

	/**
	 * @return the scoreLogrankTest
	 */
	public Double getScoreLogrankTest() {
		return scoreLogrankTest;
	}

	/**
	 * @param scoreLogrankTest the scoreLogrankTest to set
	 */
	public void setScoreLogrankTest(Double scoreLogrankTest) {
		this.scoreLogrankTest = scoreLogrankTest;
	}

	/**
	 * @return the scoreLogrankTestpvalue
	 */
	public Double getScoreLogrankTestpvalue() {
		return scoreLogrankTestpvalue;
	}

	/**
	 * @param scoreLogrankTestpvalue the scoreLogrankTestpvalue to set
	 */
	public void setScoreLogrankTestpvalue(Double scoreLogrankTestpvalue) {
		this.scoreLogrankTestpvalue = scoreLogrankTestpvalue;
	}

	/**
	 * @return the metaDataFilter
	 */
	public LinkedHashMap<String, String> getMetaDataFilter() {
		return metaDataFilter;
	}

	/**
	 * @param metaDataFilter the metaDataFilter to set
	 */
	public void setMetaDataFilter(LinkedHashMap<String, String> metaDataFilter) {
		this.metaDataFilter = metaDataFilter;
	}

	/**
	 * @return the coefficientsList
	 */
	public LinkedHashMap<String, CoxCoefficient> getCoefficientsList() {
		return coefficientsList;
	}

	/**
	 * @return the waldTestInfo
	 */
	public WaldTestInfo getWaldTestInfo() {
		return waldTestInfo;
	}

	/**
	 * @return the imat
	 */
	public double[][] getImat() {
		return imat;
	}

	/**
	 * @return the naive_imat
	 */
	public double[][] getNaive_imat() {
		return naive_imat;
	}

	/**
	 * @param degreeFreedom the degreeFreedom to set
	 */
	public void setDegreeFreedom(double degreeFreedom) {
		this.degreeFreedom = degreeFreedom;
	}

	/**
	 * @param waldTestInfo the waldTestInfo to set
	 */
	public void setWaldTestInfo(WaldTestInfo waldTestInfo) {
		this.waldTestInfo = waldTestInfo;
	}
}
