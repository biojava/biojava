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

import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Data class to represent a single sample where time and event/censor status is required
 * Additionally each variable and data associated with that variable.
 * The code handles figuring out if a variables is continuous or categorical. If categorical will
 * convert to numerical values.
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SurvivalInfo implements Comparable<SurvivalInfo> {

	private String id = "";
	private double time;
	private int status;
	private int order = 0; //not really used but included to keep track of original position if sorting.
	private double offset = 0; //offsets for linear predictor ????
	private double weight = 1; //used to set weight of survivor for over sampling.
	private int strata = 0; // this should be a boolean but leaving as an int
	private double score = 0.0;
	private double linearPredictor = 0.0;
	private double residual = 0.0;
	private String clusterValue = "";

	LinkedHashMap<String,Double> residualVariableMap = new LinkedHashMap<String,Double>();

	LinkedHashMap<String, Double> data = new LinkedHashMap<String, Double>();
	//    LinkedHashMap<String, Double> discreteData = new LinkedHashMap<String, Double>();
	LinkedHashMap<String, String> unknownDataType = new LinkedHashMap<String, String>();
	LinkedHashMap<String, String> originalMetaData = new LinkedHashMap<String,String>();

	/**
	 *
	 * @param t
	 * @param e
	 */
	public SurvivalInfo(double t, int e) {
		time = t;
		status = e;

	}

	/**
	 *
	 * @param t
	 * @param e
	 * @param d
	 */
	public SurvivalInfo(double t, int e, LinkedHashMap<String, Double> d) {
		time = t;
		status = e;

		data = d;
		for(String key : d.keySet()){
			Double value = d.get(key);
			originalMetaData.put(key, value + "");
		}
	}

	/**
	 *
	 * @param t
	 * @param e
	 * @param variable
	 * @param d
	 */
	public SurvivalInfo(double t, int e, String variable, double d) {
		time = t;
		status = e;

		data.put(variable, d);
		originalMetaData.put(variable, String.valueOf(d));
	}


	/**
	 * Set the residual value for the variable for this sample. Called from CoxScore.java
	 * @param variable
	 * @param value
	 */
	public void setResidualVariable(String variable, Double value){
		residualVariableMap.put(variable, value);
	}

	/**
	 *
	 * @param variable
	 * @return
	 */
	public Double getResidualVariable(String variable){
		return residualVariableMap.get(variable);

	}

	/**
	 *
	 * @param variable
	 * @return
	 */
	public String getUnknownDataTypeVariable(String variable){
		return unknownDataType.get(variable);
	}

	/**
	 *
	 * @param variable
	 * @return
	 */
	public String getOriginalMetaData(String variable){
		return originalMetaData.get(variable);
	}

	/**
	 *
	 * @param variable
	 * @param value
	 */
	public void addUnknownDataTypeVariable(String variable, String value) {
		originalMetaData.put(variable, value);
		unknownDataType.put(variable, value);
	}

	/**
	 *
	 * @param variable
	 * @param value
	 */
	public void updateContinousVariable(String variable, Double value){
		data.put(variable, value);
	}

	/**
	 *
	 * @param variable
	 * @param value
	 */
	public void addContinuousVariable(String variable, Double value) {
		originalMetaData.put(variable, value + "");
		data.put(variable, value);
	}

	/**
	 *
	 * @param variable
	 * @return
	 */
	public Double getContinuousVariable(String variable) {
		return data.get(variable);
	}

	/**
	 *
	 * @param groupName
	 * @return
	 */
	public ArrayList<String> getGroupCategories(String groupName) {
		ArrayList<String> groupNameList = new ArrayList<String>();
		for (String key : data.keySet()) {
			if (key.startsWith(groupName + "_")) {
				groupNameList.add(key);
			}
		}
		return groupNameList;
	}

//    public void addDiscreteVariable(String variable, double value) {
//        discreteData.put(variable, value);
//    }

//    public Double getDiscreteVariable(String variable) {
//        return discreteData.get(variable);
//    }

	/**
	 *
	 * @return
	 */
	public ArrayList<String> getDataVariables(){
		ArrayList<String> v = new ArrayList<String>();
		v.addAll(data.keySet());
		v.addAll(unknownDataType.keySet());

		return v;
	}

	/**
	 *
	 * @return
	 */
	public int getNumberVariables(){
		return data.size();
	}

	/**
	 *
	 * @param variable
	 * @return
	 */
	public Double getVariable(String variable) {
		Double value = data.get(variable);

		return value;
	}

	@Override
	public String toString() {
		return "t=" + time + " e=" + status + " o=" + order;
	}
	//    double CompNum4Sort(double[] a, double[] b) {
	//(time - time - (status -status) /1024)
	//    return (a[0] - b[0] - (a[1] - b[1]) / 1024);
	// }

	@Override
	public int compareTo(SurvivalInfo o) {
		//    return (int) (this.time - o.time - (this.status - o.status) / 1024);
		if (time < o.time) {
			return -1;
		} else if (time > o.time) {
			return 1;
		} else {
			if (this.status == o.status) {
				return 0;
			} else if (status == 1) {
				return -1;
			} else {
				return 1;
			}
		}

	}

	/**
	 * @return the offset
	 */
	public double getOffset() {
		return offset;
	}

	/**
	 * @param offset the offset to set
	 */
	public void setOffset(double offset) {
		this.offset = offset;
	}

	/**
	 * @return the weight
	 */
	public double getWeight() {
		return weight;
	}

	/**
	 * @param weight the weight to set
	 */
	public void setWeight(double weight) {
		this.weight = weight;
	}

	/**
	 * @return the strata
	 */
	public int getStrata() {
		return strata;
	}

	/**
	 * @param strata the strata to set
	 */
	public void setStrata(int strata) {
		this.strata = strata;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * @return the linearPredictor
	 */
	public double getLinearPredictor() {
		return linearPredictor;
	}

	/**
	 * @param linearPredictor the linearPredictor to set
	 */
	public void setLinearPredictor(double linearPredictor) {
		this.linearPredictor = linearPredictor;
	}

	/**
	 * @return the residual
	 */
	public double getResidual() {
		return residual;
	}

	/**
	 * @param residual the residual to set
	 */
	public void setResidual(double residual) {
		this.residual = residual;
	}

	/**
	 * @return the clusterValue
	 */
	public String getClusterValue() {
		return clusterValue;
	}

	/**
	 * @param clusterValue the clusterValue to set
	 */
	public void setClusterValue(String clusterValue) {
		this.clusterValue = clusterValue;
	}

	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}

	/**
	 * @param id the id to set
	 */
	public void setId(String id) {
		this.id = id;
	}

	/**
	 * @return the order
	 */
	public int getOrder() {
		return order;
	}

	/**
	 * @param order the order to set
	 */
	public void setOrder(int order) {
		this.order = order;
	}

	/**
	 * @return the time
	 */
	public double getTime() {
		return time;
	}

	/**
	 * @param time the time to set
	 */
	public void setTime(double time) {
		this.time = time;
	}

	/**
	 * @return the status
	 */
	public int getStatus() {
		return status;
	}

	/**
	 * @param status the status to set
	 */
	public void setStatus(int status) {
		this.status = status;
	}
}
