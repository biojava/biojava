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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 * Used to work with SurvivalInfo
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SurvivalInfoHelper {

	/**
	 * For each analysis this allows outputing of the data used in the calculations to a printstream/file. This then
	 * allows the file to be loaded into R and calculations can be verified.
	 * @param DataT
	 * @param ps
	 * @param delimiter
	 */
	public static void dump(ArrayList<SurvivalInfo> DataT, PrintStream ps, String delimiter) {
		ArrayList<String> variables = DataT.get(0).getDataVariables();
		ps.print("Seq" + delimiter);
		for (String variable : variables) {
			ps.print(variable + delimiter);
		}
		ps.print("TIME" + delimiter + "STATUS" + delimiter + "WEIGHT" + delimiter + "STRATA");

		ps.println();
		for (SurvivalInfo si : DataT) {
			ps.print(si.getOrder() + delimiter);
			for (String variable : variables) {
				Double value = si.getVariable(variable);
				ps.print(value + delimiter);
			}

			ps.print(si.getTime() + delimiter + si.getStatus() + delimiter + si.getWeight() + delimiter + si.getStrata());

			ps.println();
		}


	}

	/**
	 * If any not numeric value then categorical
	 * @param values
	 * @return
	 */
	private static boolean isCategorical(LinkedHashMap<String, Double> values) {
		try {
			for (String value : values.keySet()) {
				Double.parseDouble(value);
			}
			return false;
		} catch (Exception e) {
			return true;
		}

	}

	/**
	 * Take a collection of categorical data and convert it to numeric to be used in cox calculations
	 * @param DataT
	 */
	public static void categorizeData(ArrayList<SurvivalInfo> DataT) {

		//Go through and get all variable value pairs
		LinkedHashMap<String, LinkedHashMap<String, Double>> valueMap = new LinkedHashMap<String, LinkedHashMap<String, Double>>();
		for (SurvivalInfo si : DataT) {

			for (String key : si.unknownDataType.keySet()) {
				LinkedHashMap<String, Double> map = valueMap.get(key);
				if (map == null) {
					map = new LinkedHashMap<String, Double>();
					valueMap.put(key, map);
				}
				map.put(si.unknownDataType.get(key), null);
			}
		}

		for (String variable : valueMap.keySet()) {
			LinkedHashMap<String, Double> values = valueMap.get(variable);
			if (isCategorical(values)) {
				ArrayList<String> categories = new ArrayList<String>(values.keySet());
				Collections.sort(categories); //go ahead and put in alphabetical order
				if (categories.size() == 2) {
					for (String value : values.keySet()) {
						int index = categories.indexOf(value);
						values.put(value, index + 0.0);
					}
				} else {
					for (String value : values.keySet()) {
						int index = categories.indexOf(value);
						values.put(value, index + 1.0);
					}
				}

			} else {
				for (String value : values.keySet()) {
					Double d = Double.parseDouble(value);
					values.put(value, d);
				}
			}
		}

		for (SurvivalInfo si : DataT) {
			for (String key : si.unknownDataType.keySet()) {
				LinkedHashMap<String, Double> map = valueMap.get(key);
				String value = si.unknownDataType.get(key);
				Double d = map.get(value);
				si.data.put(key, d);
			}
		}

		for (SurvivalInfo si : DataT) {
			si.unknownDataType.clear();
		}

	}

	/**
	 * To test for interactions use two variables and create a third variable where the two are multiplied together.
	 * @param variable1
	 * @param variable2
	 * @param survivalInfoList
	 * @return
	 */
	public static ArrayList<String> addInteraction(String variable1, String variable2, ArrayList<SurvivalInfo> survivalInfoList) {
		ArrayList<String> variables = new ArrayList<String>();
		variables.add(variable1);
		variables.add(variable2);
		variables.add(variable1 + ":" + variable2);
		for (SurvivalInfo si : survivalInfoList) {
			Double value1 = si.getVariable(variable1);
			Double value2 = si.getVariable(variable2);
			Double value3 = value1 * value2;
			si.addContinuousVariable(variable1 + ":" + variable2, value3);
		}
		return variables;
	}

	/**
	 * Need to allow a range of values similar to cut in R and a continuous c
	 *
	 * @param range
	 * @param variable
	 * @param groupName
	 * @param survivalInfoList
	 * @throws Exception
	 */
	public static void groupByRange(double[] range, String variable, String groupName, ArrayList<SurvivalInfo> survivalInfoList) throws Exception {
		ArrayList<String> labels = new ArrayList<String>();
		for (int i = 0; i < range.length; i++) {
			String label = "";
			if (i == 0) {
				label = "[<=" + range[i] + "]";
			} else if (i == range.length - 1) {
				label = "[" + (range[i - 1] + 1) + "-" + range[i] + "]";
				labels.add(label);
				label = "[>" + range[i] + "]";
			} else {
				label = "[" + (range[i - 1] + 1) + "-" + range[i] + "]";
			}
			labels.add(label);
		}
		ArrayList<String> validLabels = new ArrayList<String>();

		//need to find the categories so we can set 1 and 0 and not include ranges with no values
		for (SurvivalInfo si : survivalInfoList) {
			Double value = si.getContinuousVariable(variable);
			if (value == null) {
				throw new Exception("Variable " + variable + " not found in " + si.toString());
			}
			int rangeIndex = getRangeIndex(range, value);
			String label = labels.get(rangeIndex);
			if (!validLabels.contains(groupName + "_" + label)) {
				validLabels.add(groupName + "_" + label);
			}
		}
		Collections.sort(validLabels);
		System.out.println("Valid Lables:" + validLabels);
		for (SurvivalInfo si : survivalInfoList) {
			Double value = si.getContinuousVariable(variable);
			if (value == null) {
				throw new Exception("Variable " + variable + " not found in " + si.toString());
			}
			int rangeIndex = getRangeIndex(range, value);
			String label = labels.get(rangeIndex);
			String inLable = groupName + "_" + label;
			for (String gl : validLabels) {
				if (gl.equals(inLable)) {
					si.addContinuousVariable(gl, 1.0);
				} else {
					si.addContinuousVariable(gl, 0.0);
				}
			}
		}

	}

	/**
	 *
	 * @param groupName
	 * @param survivalInfoList
	 * @return
	 */
	public static ArrayList<String> getGroupCategories(String groupName, ArrayList<SurvivalInfo> survivalInfoList) {
		return survivalInfoList.get(0).getGroupCategories(groupName);
	}

	private static int getRangeIndex(double[] range, double value) throws Exception {
		for (int i = 0; i < range.length; i++) {
			if (i == 0 && value <= range[i]) {
				return i;
			}
			if (value <= range[i]) {
				return i;
			}

		}

		if (value > range[range.length - 1]) {
			return range.length;
		}
		throw new Exception("Value " + value + " not found in range ");
	}
}
