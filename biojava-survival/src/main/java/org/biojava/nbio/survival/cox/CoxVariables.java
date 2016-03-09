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

import java.text.DecimalFormat;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxVariables {

	//   public GeneSetResults gsr;
	private String cohortName = "";
	private String geneSet = "";
	private String genes = "";

	/**
	 *
	 * @param cohortName
	 * @param geneSet
	 * @param genes
	 */
	public CoxVariables(String cohortName, String geneSet, String genes) {
		this.cohortName = cohortName;
		this.geneSet = geneSet;
		this.genes = genes;
	}

	/**
	 * Need a unique id from String
	 *
	 * @return
	 */
	public int getUniqueID() {
		String link = geneSet + "_" + cohortName;
		return link.hashCode();
	}
	private LinkedHashMap<String, CoxInfo> coxInfoHashMap = new LinkedHashMap<String, CoxInfo>();

	/**
	 *
	 * @param name
	 * @param coxInfo
	 */
	public void putCoxInfo(String name, CoxInfo coxInfo) {
		coxInfoHashMap.put(name, coxInfo);

	}

	/**
	 *
	 * @param name
	 * @return
	 */
	public CoxInfo getCoxInfo(String name) {
		return coxInfoHashMap.get(name);

	}

	/**
	 *
	 * @param file
	 * @return
	 */
	public String encodeFileURL(String file) {
		file = file.replaceAll(" ", "%20");
		file = file.replaceAll("<", "%3C");
		file = file.replaceAll(">", "%3E");
		return file;
	}
	//   static GeneProfiler geneProfiler = null;

	/**
	 *
	 * @param d
	 * @param precision
	 * @param pad
	 * @return
	 */
	public static String fmt(Double d, int precision, int pad) {
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

	@Override
	public String toString() {
		String coxOutput = geneSet + "\r\n";
		//    co = co + genes + "\r\n";
		coxOutput = coxOutput + cohortName + "," + genes.replace(',', ' ') + "\r\n";
		coxOutput = coxOutput + ",Coe,StdErr,p-value,HR,HR Lo 95%,HR Hi 95%\r\n";
		for (String variables : coxInfoHashMap.keySet()) {
			CoxInfo ci = coxInfoHashMap.get(variables);

			coxOutput = coxOutput + "Overall Model Fit p-value=" + fmt(ci.getOverallModelFitPvalue(), 5, 0) + "\r\n";
			coxOutput = coxOutput + ci.getCoefficientText(false, "", ",", "", "");
			coxOutput = coxOutput + "\r\n";


		}

		return coxOutput;
	}

	/**
	 * @return the cohortName
	 */
	public String getCohortName() {
		return cohortName;
	}

	/**
	 * @return the geneSet
	 */
	public String getGeneSet() {
		return geneSet;
	}

	/**
	 * @return the genes
	 */
	public String getGenes() {
		return genes;
	}

	/**
	 * @return the coxInfoHashMap
	 */
	public LinkedHashMap<String, CoxInfo> getCoxInfoHashMap() {
		return coxInfoHashMap;
	}
}
