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
package org.biojava.nbio.survival.kaplanmeier.figure;

import org.biojava.nbio.survival.cox.StrataInfo;
import org.biojava.nbio.survival.cox.SurvFitInfo;
import org.biojava.nbio.survival.cox.SurvivalInfo;
import org.biojava.nbio.survival.data.WorkSheet;

import javax.swing.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 * Ported from survfitKM.S When combining multiple entries with same time not
 * sure how the weighting adds up
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SurvFitKM {

	/**
	 *
	 */
	public enum Method {

		/**
		 *
		 */
		kaplanMeier,
		/**
		 *
		 */
		flemingHarrington,
		/**
		 *
		 */
		fh2;
	}

	/**
	 *
	 */
	public enum Error {

		/**
		 *
		 */
		greenwood,
		/**
		 *
		 */
		tsiatis;
	}

	/**
	 *
	 */
	public enum ConfType {

		/**
		 *
		 */
		log,
		/**
		 *
		 */
		log_log,
		/**
		 *
		 */
		plain,
		/**
		 *
		 */
		none;
	}

	/**
	 *
	 */
	public enum ConfLower {

		/**
		 *
		 */
		usual,
		/**
		 *
		 */
		peto,
		/**
		 *
		 */
		modified;
	}

	/**
	 *
	 * @param survivalData
	 * @param useWeights
	 * @return
	 * @throws Exception
	 */
	public SurvFitInfo process(LinkedHashMap<String, ArrayList<CensorStatus>> survivalData, boolean useWeights) throws Exception {
		ArrayList<SurvivalInfo> survivalInfoList = new ArrayList<SurvivalInfo>();
		int i = 0;
		for (String strata : survivalData.keySet()) {
			ArrayList<CensorStatus> csList = survivalData.get(strata);
			for (CensorStatus cs : csList) {
				SurvivalInfo si = new SurvivalInfo(cs.time, Integer.parseInt(cs.censored));
				si.setOrder(i);
				i++;
				si.setWeight(cs.weight);
				si.addUnknownDataTypeVariable("STRATA", strata);
				si.addUnknownDataTypeVariable("VALUE", cs.value + "");
				survivalInfoList.add(si);
			}
		}

		return process("STRATA", survivalInfoList, useWeights);

	}

	/**
	 *
	 * @param datafile
	 * @param timeColumn
	 * @param statusColumn
	 * @param weightColumn
	 * @param variableColumn
	 * @param useWeights
	 * @return
	 * @throws Exception
	 */
	public SurvFitInfo process(String datafile, String timeColumn, String statusColumn, String weightColumn, String variableColumn, boolean useWeights) throws Exception {
		WorkSheet worksheet = WorkSheet.readCSV(datafile, '\t');
		ArrayList<SurvivalInfo> survivalInfoList = new ArrayList<SurvivalInfo>();

		int i = 1;
		for (String row : worksheet.getRows()) {
			double time = worksheet.getCellDouble(row, timeColumn);

			double c = worksheet.getCellDouble(row, statusColumn);
			double weight = 1.0;
			if (weightColumn != null && weightColumn.length() > 0) {
				weight = worksheet.getCellDouble(row, weightColumn);

			}
			int strata = 0;

			int censor = (int) c;

			if (weight <= 0) {
				//   System.out.println("Weight <= 0 Sample=" + row + " weight=" + weight);
				i++;
				continue;
			}



			SurvivalInfo si = new SurvivalInfo(time, censor);
			si.setOrder(i);
			si.setWeight(weight);
			si.setStrata(strata);


			String value = worksheet.getCell(row, variableColumn);
			si.addUnknownDataTypeVariable(variableColumn, value);



			survivalInfoList.add(si);
			i++;
		}



		return process(variableColumn, survivalInfoList, useWeights);
	}

	/**
	 *
	 * @param variable
	 * @param dataT
	 * @param useWeighted
	 * @return
	 * @throws Exception
	 */
	public SurvFitInfo process(String variable, ArrayList<SurvivalInfo> dataT, boolean useWeighted) throws Exception {
		return this.process(variable, dataT, Method.kaplanMeier, Error.greenwood, true, .95, ConfType.log, ConfLower.usual, null, null, useWeighted);
	}


	public LinkedHashMap<String, StrataInfo>  processStrataInfo(String variable, ArrayList<SurvivalInfo> dataT, SurvFitKM.Method method, SurvFitKM.Error error, boolean seFit, double confInt, ConfType confType, ConfLower confLower, Double startTime, Double newTime, boolean useWeighted) throws Exception{
				Collections.sort(dataT);
		if (startTime == null && newTime != null) {
			startTime = newTime;
		}
		//int ny = 2; // setup for right censored versus counting
		if (startTime != null) {
			throw new Exception("Filter on startTime not implemented");
		}
		int n = dataT.size();


		LinkedHashMap<String, Integer> levels = new LinkedHashMap<String, Integer>();
		LinkedHashMap<String, ArrayList<SurvivalInfo>> strataHashMap = new LinkedHashMap<String, ArrayList<SurvivalInfo>>();

		for (int i = 0; i < n; i++) {
			SurvivalInfo si = dataT.get(i);
			String value = si.getUnknownDataTypeVariable(variable);
			Integer count = levels.get(value);
			if (count == null) {
				count = 0;
			}
			count++;
			levels.put(value, count);
			ArrayList<SurvivalInfo> strataList = strataHashMap.get(value);
			if (strataList == null) {
				strataList = new ArrayList<SurvivalInfo>();
				strataHashMap.put(value, strataList);
			}
			strataList.add(si);

		}

		//int nstrat = levels.size();

		LinkedHashMap<String, StrataInfo> strataInfoHashMap = new LinkedHashMap<String, StrataInfo>();

		for (String strata : strataHashMap.keySet()) {

			ArrayList<SurvivalInfo> strataList = strataHashMap.get(strata);
			StrataInfo strataInfo = new StrataInfo();
			strataInfoHashMap.put(strata, strataInfo);


			Double previousTime = null;
			for (SurvivalInfo si : strataList) {
				double w = 1.0;
				if (useWeighted) {
					w = si.getWeight();
				}

				if (previousTime == null || si.getTime() != previousTime) {
					strataInfo.getTime().add(si.getTime());
					if (si.getStatus() == 0) {
						strataInfo.getStatus().add(0);
						strataInfo.getNcens().add(w);
						strataInfo.getNevent().add(0.0);
					} else {
						strataInfo.getNcens().add(0.0);
						strataInfo.getNevent().add(w);
						strataInfo.getStatus().add(1);
					}
					strataInfo.getNrisk().add(0.0);
					strataInfo.getStderr().add(0.0);
					strataInfo.getWeight().add(w);
				} else {
					//we have the same time so add to previous entry
					int index = strataInfo.getTime().size() - 1;
					if (si.getStatus() == 0) {
						double nw = strataInfo.getNcens().get(index) + w;
						strataInfo.getNcens().remove(index);
						strataInfo.getNcens().add(nw);
						// strataInfo.nevent.add(0.0);
					} else {
						// strataInfo.ncens.add(0.0);
						double nw = strataInfo.getNevent().get(index) + w;
						strataInfo.getNevent().remove(index);
						strataInfo.getNevent().add(nw);

					}
					double nw = strataInfo.getWeight().get(index) + w;
					strataInfo.getWeight().remove(index);
					strataInfo.getWeight().add(nw);
				}
				previousTime = si.getTime();
				//  strataInfo.status.add(si.status);



				Integer ndead = strataInfo.getNdead().get(si.getTime());
				if (ndead == null) {
					ndead = 0;
				}
				if (si.getStatus() == 1) {
					ndead++;
				}
				strataInfo.getNdead().put(si.getTime(), ndead);

			}



			int j = strataInfo.getWeight().size() - 1;
			double cw = 0.0;
			for (int i = strataInfo.getWeight().size() - 1; i >= 0; i--) {
				double c = strataInfo.getWeight().get(i);
				cw = cw + c;
				strataInfo.getNrisk().set(j, cw);
				j--;
			}
			if (method == Method.kaplanMeier) {

				for (int i = 0; i < strataInfo.getNrisk().size(); i++) {
					double t = (strataInfo.getNrisk().get(i) - strataInfo.getNevent().get(i)) / strataInfo.getNrisk().get(i);
					if (i == 0) {
						strataInfo.getSurv().add(t);
					} else {
						strataInfo.getSurv().add(t * strataInfo.getSurv().get(i - 1));
					}
				}
			} else if (method == Method.flemingHarrington) {
				for (int i = 0; i < strataInfo.getNrisk().size(); i++) {
					double hazard = (strataInfo.getNevent().get(i)) / strataInfo.getNrisk().get(i);
					if (i == 0) {
						strataInfo.getSurv().add(Math.exp(-1.0 * hazard));
					} else {
						strataInfo.getSurv().add(Math.exp(-1.0 * (hazard + strataInfo.getSurv().get(i - 1))));
					}
				}
			} else if (method == Method.fh2) {
				throw new Exception("Method.fh2 not supported. Need to implement survfit4.c ");
			}

			if (seFit) {
				if (error == Error.greenwood) {
					for (int i = 0; i < strataInfo.getNrisk().size(); i++) {
						double t = (strataInfo.getNevent().get(i)) / (strataInfo.getNrisk().get(i) * (strataInfo.getNrisk().get(i) - strataInfo.getNevent().get(i)));
						if (i == 0) {
							strataInfo.getVarhaz().add(t);
						} else {
							strataInfo.getVarhaz().add(t + strataInfo.getVarhaz().get(i - 1));
						}
					}
				} else if (method == Method.kaplanMeier || method == Method.flemingHarrington) {
					for (int i = 0; i < strataInfo.getNrisk().size(); i++) {
						double t = (strataInfo.getNevent().get(i)) / (strataInfo.getNrisk().get(i) * strataInfo.getNrisk().get(i));
						if (i == 0) {
							strataInfo.getVarhaz().add(t);
						} else {
							strataInfo.getVarhaz().add(t + strataInfo.getVarhaz().get(i - 1));
						}
					}

				} else {
					//varhaz[[i]] <- cumsum(nevent* tsum$sum2)
					throw new Exception("Method.fh2 not supported. Need to implement survfit4.c ");
				}
				strataInfo.getStderr().clear();
				for (int i = 0; i < strataInfo.getNrisk().size(); i++) {

					strataInfo.getStderr().add(Math.sqrt(strataInfo.getVarhaz().get(i)));
				}


			}

		}

		ArrayList<Boolean> events = new ArrayList<Boolean>();
		ArrayList<Double> nrisk = new ArrayList<Double>();
		for (StrataInfo strataInfo : strataInfoHashMap.values()) {
			boolean firsttime = true;
			for (int j = 0; j < strataInfo.getNevent().size(); j++) {
				Double d = strataInfo.getNevent().get(j);
				if (d > 0 || firsttime) {
					firsttime = false;
					events.add(true);
					nrisk.add(strataInfo.getNrisk().get(j));
				} else {
					events.add(false);
				}
			}

		}
		ArrayList<Integer> zz = new ArrayList<Integer>();
		for (int i = 0; i < events.size(); i++) {
			if (events.get(i)) {
				zz.add(i + 1);
			}
		}
		zz.add(events.size() + 1);
		ArrayList<Integer> diffzz = new ArrayList<Integer>();
		for (int i = 0; i < zz.size() - 1; i++) {
			diffzz.add(zz.get(i + 1) - zz.get(i));
		}
		//System.out.println(diffzz);
		ArrayList<Double> nlag = new ArrayList<Double>();
		for (int j = 0; j < nrisk.size(); j++) {
			int count = diffzz.get(j);
			for (int c = 0; c < count; c++) {
				nlag.add(nrisk.get(j));
			}
		}

		// System.out.println(nlag);
		// System.out.println("nlag.size=" + nlag.size());
		if (confLower == ConfLower.usual) {
			for (StrataInfo strataInfo : strataInfoHashMap.values()) {
				strataInfo.setStdlow(strataInfo.getStderr());
			}
		} else if (confLower == ConfLower.peto) {
			for (StrataInfo strataInfo : strataInfoHashMap.values()) {
				for (int j = 0; j < strataInfo.getSurv().size(); j++) {
					double v = Math.sqrt((1.0 - strataInfo.getSurv().get(j)) / strataInfo.getNrisk().get(j));
					strataInfo.getStdlow().add(v);
				}
			}
		} else if (confLower == ConfLower.modified) {
			int i = 0;
			for (StrataInfo strataInfo : strataInfoHashMap.values()) {
				for (int j = 0; j < strataInfo.getSurv().size(); j++) {
					double v = strataInfo.getStderr().get(j) * Math.sqrt(nlag.get(i) / strataInfo.getNrisk().get(j));
					strataInfo.getStdlow().add(v);
					i++;
				}
			}
		}

		//zval <- qnorm(1- (1-conf.int)/2, 0,1)
		double zvalue = 1.959964;

		if (confType == ConfType.plain) {
			for (StrataInfo strataInfo : strataInfoHashMap.values()) {
				for (int j = 0; j < strataInfo.getSurv().size(); j++) {
					double upper = strataInfo.getSurv().get(j) + zvalue * strataInfo.getStderr().get(j) * strataInfo.getSurv().get(j);
					double lower = strataInfo.getSurv().get(j) - zvalue * strataInfo.getStdlow().get(j) * strataInfo.getSurv().get(j);
					strataInfo.getLower().add(lower);
					strataInfo.getUpper().add(upper);
				}
			}

		} else if (confType == ConfType.log) {
			for (StrataInfo strataInfo : strataInfoHashMap.values()) {
				for (int j = 0; j < strataInfo.getSurv().size(); j++) {
					double surv = strataInfo.getSurv().get(j);
					if (surv == 0) {
						strataInfo.getLower().add(Double.NaN);
						strataInfo.getUpper().add(Double.NaN);
					} else {
						double upper = Math.exp(Math.log(surv) + zvalue * strataInfo.getStderr().get(j));
						double lower = Math.exp(Math.log(surv) - zvalue * strataInfo.getStdlow().get(j));
						strataInfo.getLower().add(lower);
						strataInfo.getUpper().add(upper);
					}
				}
			}

		} else if (confType == ConfType.log_log) {
			throw new Exception("ConfType log-log currently not supported");
		}


//        if (false) {
//            for (String strata : strataInfoHashMap.keySet()) {
//                StrataInfo strataInfo = strataInfoHashMap.get(strata);
//                System.out.println(strataInfo.toString());
//                System.out.println();
//            }
//            System.out.println();
//        }

		return strataInfoHashMap;
	}

	/**
	 *
	 * @param variable
	 * @param dataT
	 * @param method
	 * @param error
	 * @param seFit
	 * @param confInt
	 * @param confType
	 * @param confLower
	 * @param startTime
	 * @param newTime
	 * @param useWeighted
	 * @return
	 * @throws Exception
	 */
	public SurvFitInfo process(String variable, ArrayList<SurvivalInfo> dataT, SurvFitKM.Method method, SurvFitKM.Error error, boolean seFit, double confInt, ConfType confType, ConfLower confLower, Double startTime, Double newTime, boolean useWeighted) throws Exception {
		SurvFitInfo si = new SurvFitInfo();

		LinkedHashMap<String, StrataInfo> strataInfoHashMap = this.processStrataInfo(variable, dataT, method, error, seFit, confInt, confType, confLower, startTime, newTime, useWeighted);
		si.setStrataInfoHashMap(strataInfoHashMap);
		LinkedHashMap<String, StrataInfo> unweightedStrataInfoHashMap = this.processStrataInfo(variable, dataT, method, error, seFit, confInt, confType, confLower, startTime, newTime, false);
		si.setUnweightedStrataInfoHashMap(unweightedStrataInfoHashMap);
		si.setWeighted(useWeighted);
		return si;
	}

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
		try {
			String datafile = "/Users/Scooter/scripps/ngs/BLJ/E2197/Predictive Signatures/V$HSF_Q6-E2197 TTR.txt";

			SurvFitKM survFitKM = new SurvFitKM();

			SurvFitInfo si = survFitKM.process(datafile, "TIME", "STATUS", "WEIGHT", "MEAN", true);



			if (true) {

				KaplanMeierFigure kaplanMeierFigure = new KaplanMeierFigure();





				JFrame application = new JFrame();
				application.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				application.add(kaplanMeierFigure);
				kaplanMeierFigure.setSize(500, 400);

				application.setSize(500, 400);         // window is 500 pixels wide, 400 high
				application.setVisible(true);

				ArrayList<String> titles = new ArrayList<String>();
				titles.add("Line 1");
				titles.add("line 2");
				kaplanMeierFigure.setSurvivalData(titles, si, null);

				ArrayList<String> figureInfo = new ArrayList<String>();
				//   figureInfo.add("HR=2.1 95% CI(1.8-2.5)");
				//   figureInfo.add("p-value=.001");
				kaplanMeierFigure.setFigureLineInfo(figureInfo);

				kaplanMeierFigure.savePNG("/Users/Scooter/Downloads/test.png");



			}



		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
