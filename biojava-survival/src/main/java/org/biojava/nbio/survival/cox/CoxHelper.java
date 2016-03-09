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


import org.biojava.nbio.survival.data.WorkSheet;

import java.util.ArrayList;

/**
 * The CoxHelper class is provided to start with a tab delimited file in a similar process in R and return the results as a CoxInfo class.
 * Given the number of options for adjusting the calculations using weighting, strata, clustering etc the helper class can be used to hide
 * the complexity for typical use case.
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxHelper {

	/**
	 *
	 * @param datafile The tab delimited file containing survival data and variables. The first column needs to be unique index
	 * @param timeColumn The column representing the event/censor time
	 * @param statusColumn The column representing an event=1 and censor=0
	 * @param weightColumn For case-cohort data sets may require weighting to reflect the entire cohort
	 * @param strataColumn A column representing strata data
	 * @param clusterColumn If robost variation calculation is required the cluster column will group samples by the value in this column
	 * @param variables The variables to be used in the cox regression analysis. For Interactions using variable1:variable2
	 * @param useStrata Boolean to indicate if strata column should be used
	 * @param useWeights Boolean to indicate if weight column should be used
	 * @return
	 * @throws Exception
	 */


	public static CoxInfo process(String datafile, String timeColumn, String statusColumn, String weightColumn, String strataColumn, String clusterColumn, ArrayList<String> variables, boolean useStrata, boolean useWeights) throws Exception {
		WorkSheet worksheet = WorkSheet.readCSV(datafile, '\t');
		return process(worksheet, timeColumn, statusColumn, weightColumn, strataColumn, clusterColumn, variables, useStrata, useWeights);
	}

	/**
	 *
	 * @param worksheet
	 * @param timeColumn The column representing the event/censor time
	 * @param statusColumn The column representing an event=1 and censor=0
	 * @param weightColumn For case-cohort data sets may require weighting to reflect the entire cohort
	 * @param strataColumn A column representing strata data
	 * @param clusterColumn If robost variation calculation is required the cluster column will group samples by the value in this column
	 * @param variables The variables to be used in the cox regression analysis. For Interactions using variable1:variable2
	 * @param useStrata Boolean to indicate if strata column should be used
	 * @param useWeights Boolean to indicate if weight column should be used
	 * @return
	 */
	public static CoxInfo process(WorkSheet worksheet, String timeColumn, String statusColumn, String weightColumn, String strataColumn, String clusterColumn, ArrayList<String> variables, boolean useStrata, boolean useWeights) {

		try {
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
				if (strataColumn != null && strataColumn.length() > 0) {
					strata = worksheet.getCellDouble(row, strataColumn).intValue();
				}
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
				for (String column : variables) {
					if (column.contains(":")) {
						continue;
					}
					String value = worksheet.getCell(row, column);
					si.addUnknownDataTypeVariable(column, value);
				}
				if (clusterColumn != null && clusterColumn.length() > 0) {
					String v = worksheet.getCell(row, clusterColumn);
					si.setClusterValue(v);
				}

				survivalInfoList.add(si);
				i++;
			}



			boolean cluster = false;
			boolean robust = false;
			if (clusterColumn != null && clusterColumn.length() > 0) {
				cluster = true;
				robust = true;
			}
			//       variables.add("TREAT:AGE");
			CoxR cox = new CoxR();
			CoxInfo ci = cox.process(variables, survivalInfoList, useStrata, useWeights, robust, cluster);
			// System.out.println(ci);

			//applying Bob Gray's correction for weighted strata wtexamples.docx
			//           CoxCC.process(ci, survivalInfoList);
			//           ci.dump();
			//           ci.calcSummaryValues();

			return ci;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}



	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		// TODO code application logic here
		try {
			if (true) {
				String datafile = "/Users/Scooter/scripps/ngs/DataSets/E2197/misc/ecoglabtransfer/500790/2013.05.10.12.28.58.313/clindasl0228.txt";
				ArrayList<String> variables = new ArrayList<String>();
				variables.add("nndpos");
				variables.add("meno");
//              variables.add("er1");
//              variables.add("meno:er1");

				CoxInfo ci = CoxHelper.process(datafile, "ttr", "recind", "wt", "sstrat", "Seq", variables, false, true);

			  //  ci.dump();

				System.out.println(ci);
				System.out.println();

				CoxCC.process(ci);

				ci.dump();

			}


		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
