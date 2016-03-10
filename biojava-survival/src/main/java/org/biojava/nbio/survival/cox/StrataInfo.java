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
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Information needed to represent a survival curve
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class StrataInfo {

	private ArrayList<Double> time = new ArrayList<Double>();
	private ArrayList<Integer> status = new ArrayList<Integer>();
	private ArrayList<Double> nevent = new ArrayList<Double>();
	private ArrayList<Double> ncens = new ArrayList<Double>();
	private ArrayList<Double> nrisk = new ArrayList<Double>();
	private ArrayList<Double> weight = new ArrayList<Double>();
	private ArrayList<Double> surv = new ArrayList<Double>();
	private ArrayList<Double> varhaz = new ArrayList<Double>();
	private ArrayList<Double> stderr = new ArrayList<Double>();
	private ArrayList<Double> stdlow = new ArrayList<Double>();
	private ArrayList<Double> upper = new ArrayList<Double>();
	private ArrayList<Double> lower = new ArrayList<Double>();
	private LinkedHashMap<Double, Integer> ndead = new LinkedHashMap<Double, Integer>();
	DecimalFormat df = new DecimalFormat("#.######");
	DecimalFormat dfe = new DecimalFormat("0.000000E0");

	/**
	 * Need to find the actual time for the nearest time represented as a
	 * percentage Would be used to then look up the number at risk at that
	 * particular time
	 *
	 * @param timePercentage
	 * @return
	 */
	public Double getNearestTime(double timePercentage) {
		//the arrays should be sorted by time so this step is probably not needed
		Double minTime = null;
		Double maxTime = null;
		for (Double t : time) {
			if (minTime == null || t < minTime) {
				minTime = t;
			}
			if (maxTime == null || t > maxTime) {
				maxTime = t;
			}
		}
		Double timeRange = maxTime - minTime;
		Double targetTime = minTime + timePercentage * timeRange;
		Double previousTime = null;
		for (Double t : time) {
			if (previousTime == null || t <= targetTime) {
				previousTime = t;
			} else {
				return previousTime;
			}
		}
		return previousTime;
	}

	/**
	 * Selection of number of risk will depend on the precision and rounding of
	 * time in the survival table. If you are asking for 12 and entry exists for
	 * 11.9999999 then 12 is greater than 11.99999 unless you round.
	 *
	 * @param t
	 * @return
	 */
	public Double getNearestAtRisk(double t) {
		Integer index = 0;
/*       String timeValue = t + "";
		String format = "#";
		int numDecimals = 0;
		int decimalIndex = timeValue.indexOf(".");
		if (decimalIndex > 0) {
			for (int i = timeValue.length() - 1; i > decimalIndex; i--) {
				if (timeValue.charAt(i) == '0' && numDecimals == 0) {
					continue;
				}
				if (i == decimalIndex - 1) {
					format = format + ".#";
				} else {
					format = format + "#";
				}
			}
		}
 */
		DecimalFormat newFormat = new DecimalFormat("#.#"); //used to round on expected precision of time. Not correct but trying to match the other packages

		for (int i = 0; i < time.size(); i++) {
			Double compareTime = time.get(i);
		  //  compareTime = new Double(Math.round(compareTime)); //this is rounding up so that we stop on the first match trying to get this to match another report. Not correct or the other report is wrong
			compareTime = Double.valueOf(newFormat.format(compareTime));
			if (compareTime < t) {
				index = i + 1 ;
			} else if(compareTime == t){
				index = i;
				break;
			}else {
				break;
			}
		}

		//http://www.inside-r.org/packages/cran/rms/docs/survplot
		//per validation using survplot from RMS package and ggkm they select the next
		//time in the future which doesn't seem to be correct as the next time represents
		//knowledge about the future but maybe nrisk at that point in time is defined
		//as the nrisk prior to that time. This appears to be the case where at time 0
		//you would expect that everyone is at risk and you should report that time which
		//is the case in survplot. Added in index = 0 or if the time you are requesting has
		//an exact match
		//survplot(kma,n.risk=TRUE,time.inc=1090)
		//ggkm(kma,timeby=1090)
	   //     if(index != 0 && time.get(index) != t){
	   //      index++;
	   //     }
		if (index >= nrisk.size()) {
			return null;
		} else {
			return nrisk.get(index);
		}
	}

	/**
	 *
	 * @param d
	 * @return
	 */
	public String f(Double d) {
		String v = df.format(d);
		int l = 10 - v.length();
		for (int i = 0; i < l; i++) {
			v = v + " ";
		}
		return v;
	}

	@Override
	public String toString() {
		String o = "";
		o = o + "n=" + nevent.size() + "\r\n";
		o = o + "     time      nevent     ncens     nrisk     weight     surv   varhaz    stderr    stdlow    lower    upper\r\n";
		for (int i = 0; i < nevent.size(); i++) {
			//    if(nevent.get(i) == 0)
			//        continue;
			o = o + (i + 1) + "    " + f(time.get(i)) + " " + f(nevent.get(i)) + " " + f(ncens.get(i)) + " " + f(nrisk.get(i)) + " " + f(weight.get(i)) + " " + f(surv.get(i)) + " " + (varhaz.get(i)) + "  " + stderr.get(i) + "  " + stdlow.get(i) + "  " + lower.get(i) + "  " + upper.get(i) + "\r\n";

		}
		o = o + "\r\n";
		//   for(Integer i : ndead.values()){
		//       o = o + i + "\r\n";
		//   }

		return o;
	}

	/**
	 * @return the time
	 */
	public ArrayList<Double> getTime() {
		return time;
	}

	/**
	 * @return the surv
	 */
	public ArrayList<Double> getSurv() {
		return surv;
	}

	/**
	 * @return the stderr
	 */
	public ArrayList<Double> getStderr() {
		return stderr;
	}

	/**
	 * @return the upper
	 */
	public ArrayList<Double> getUpper() {
		return upper;
	}

	/**
	 * @return the lower
	 */
	public ArrayList<Double> getLower() {
		return lower;
	}

	/**
	 * @return the status
	 */
	public ArrayList<Integer> getStatus() {
		return status;
	}

	/**
	 * @return the nevent
	 */
	public ArrayList<Double> getNevent() {
		return nevent;
	}

	/**
	 * @return the ncens
	 */
	public ArrayList<Double> getNcens() {
		return ncens;
	}

	/**
	 * @return the nrisk
	 */
	public ArrayList<Double> getNrisk() {
		return nrisk;
	}

	/**
	 * @return the weight
	 */
	public ArrayList<Double> getWeight() {
		return weight;
	}

	/**
	 * @return the ndead
	 */
	public LinkedHashMap<Double, Integer> getNdead() {
		return ndead;
	}

	/**
	 * @return the varhaz
	 */
	public ArrayList<Double> getVarhaz() {
		return varhaz;
	}

	/**
	 * @return the stdlow
	 */
	public ArrayList<Double> getStdlow() {
		return stdlow;
	}

	/**
	 * @param stdlow the stdlow to set
	 */
	public void setStdlow(ArrayList<Double> stdlow) {
		this.stdlow = stdlow;
	}
}
