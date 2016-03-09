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

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxCoefficient {

	String name;
	double coeff; //beta
	double stdError; //se
	double robustStdError; //nse
	double z;
	double hazardRatio; //exp(beta)
	double hazardRatioLoCI;
	double hazardRatioHiCI;
	double pvalue;
	double mean;
	double standardDeviation;

	/**
	 *
	 */
	public CoxCoefficient() {
	}

	@Override
	public String toString() {
		return name + " " + coeff + " " + pvalue + " " + hazardRatio + " " + hazardRatioLoCI + " " + hazardRatioHiCI;
	}

	/**
	 *
	 * @return
	 */
	public String getHRText() {
		return fmt(hazardRatio, 2, 0) + " CI(" + fmt(hazardRatioLoCI, 2, 0) + "-" + fmt(hazardRatioHiCI, 2, 0) + ")";
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the coeff
	 */
	public double getCoeff() {
		return coeff;
	}

	/**
	 * @return the stdError
	 */
	public double getStdError() {
		return stdError;
	}

	/**
	 * @return the robustStdError
	 */
	public double getRobustStdError() {
		return robustStdError;
	}

	/**
	 * @return the z
	 */
	public double getZ() {
		return z;
	}

	/**
	 * @return the hazardRatio
	 */
	public double getHazardRatio() {
		return hazardRatio;
	}

	/**
	 * @return the hazardRatioLoCI
	 */
	public double getHazardRatioLoCI() {
		return hazardRatioLoCI;
	}

	/**
	 * @return the hazardRatioHiCI
	 */
	public double getHazardRatioHiCI() {
		return hazardRatioHiCI;
	}

	/**
	 * @return the pvalue
	 */
	public double getPvalue() {
		return pvalue;
	}

	/**
	 * @return the mean
	 */
	public double getMean() {
		return mean;
	}

	/**
	 * @return the standardDeviation
	 */
	public double getStandardDeviation() {
		return standardDeviation;
	}

	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @param coeff the coeff to set
	 */
	public void setCoeff(double coeff) {
		this.coeff = coeff;
	}

	/**
	 * @param stdError the stdError to set
	 */
	public void setStdError(double stdError) {
		this.stdError = stdError;
	}

	/**
	 * @param robustStdError the robustStdError to set
	 */
	public void setRobustStdError(double robustStdError) {
		this.robustStdError = robustStdError;
	}

	/**
	 * @param z the z to set
	 */
	public void setZ(double z) {
		this.z = z;
	}

	/**
	 * @param hazardRatio the hazardRatio to set
	 */
	public void setHazardRatio(double hazardRatio) {
		this.hazardRatio = hazardRatio;
	}

	/**
	 * @param hazardRatioLoCI the hazardRatioLoCI to set
	 */
	public void setHazardRatioLoCI(double hazardRatioLoCI) {
		this.hazardRatioLoCI = hazardRatioLoCI;
	}

	/**
	 * @param hazardRatioHiCI the hazardRatioHiCI to set
	 */
	public void setHazardRatioHiCI(double hazardRatioHiCI) {
		this.hazardRatioHiCI = hazardRatioHiCI;
	}

	/**
	 * @param pvalue the pvalue to set
	 */
	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}

	/**
	 * @param mean the mean to set
	 */
	public void setMean(double mean) {
		this.mean = mean;
	}

	/**
	 * @param standardDeviation the standardDeviation to set
	 */
	public void setStandardDeviation(double standardDeviation) {
		this.standardDeviation = standardDeviation;
	}

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
}
