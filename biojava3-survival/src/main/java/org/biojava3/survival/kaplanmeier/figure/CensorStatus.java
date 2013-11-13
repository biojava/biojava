/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.kaplanmeier.figure;

/**
 * 
 * @author willishf@gmail.com
 */
public class CensorStatus implements Comparable<CensorStatus> {

    /**
     *
     */
    public String row;
    /**
     *
     */
    public Double time;
    /**
     *
     */
    public String censored;
    /**
     *
     */
    public String group;
    /**
     *
     */
    public Double value;
    /**
     *
     */
    public Double zscore;
    /**
     *
     */
    public Double weight = 1.0; // assume default weight 1.0
    
    
    private Double percentage = null; //allow the percentage to be set externally for various weighted correction methods.
    /**
     *
     */
    public Double nevents;
    /**
     *
     */
    public Double ncens;
    /**
     *
     */
    public Double nrisk;

    /**
     *
     */
    public CensorStatus() {
    }

    /**
     *
     * @param group
     * @param time
     * @param censored
     */
    public CensorStatus(String group, Double time, String censored) {
        this.group = group;
        this.time = time;
        this.censored = censored;
    }

    /**
     *
     * @return
     */
    public CensorStatus getCopy(){
        CensorStatus cs = new CensorStatus();
        cs.row = row;
        cs.time = time;
        cs.censored = censored;
        cs.group = group;
        cs.value = value;
        cs.zscore = zscore;
        return cs;
    }
    
    public String toString() {
        return time + " " + censored + " " + group + " " + row;
    }

    public int compareTo(CensorStatus o) {
  //    System.out.println("Comparing " + this + " " + o);
        if (time == null) {
            return -1;
        }
        if (o.time == null) {
            return 1;
        }

        if (time < o.time) {
            return -1;
        } else if (time > o.time) {
            return 1;
        } else {
            if (censored.equals(o.censored)) {
                return 0;
            }
            if (censored.equals("0")) {
                return -1;
            } else {
                return 1;
            }
        }
    }

    /**
     * @return the percentage
     */
    public Double getPercentage() {
        return percentage;
    }

    /**
     * @param percentage the percentage to set
     */
    public void setPercentage(Double percentage) {
        this.percentage = percentage;
    }
}
