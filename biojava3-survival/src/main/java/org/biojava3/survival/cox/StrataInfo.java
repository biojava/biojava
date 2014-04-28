/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Information needed to represent a survival curve
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
   
   private LinkedHashMap<Double,Integer> ndead = new LinkedHashMap<Double,Integer>();

   DecimalFormat df = new DecimalFormat("#.######");
   DecimalFormat dfe = new DecimalFormat("0.000000E0");
    /**
     *
     * @param d
     * @return
     */
    public String f(Double d){
       String v = df.format(d);
       int l = 10 - v.length();
       for(int i = 0; i < l; i++){
           v = v + " ";
       }
       return v;
   }
   
    @Override
    public String toString() {
        String o = "";
        o = o + "n=" + nevent.size() + "\r\n";
        o = o + "     time      nevent     ncens     nrisk     weight     surv   varhaz    stderr    stdlow    lower    upper\r\n"; 
        for(int i = 0; i < nevent.size(); i++){
            o = o + (i+1) + "    " + f(time.get(i)) + " " + f(nevent.get(i))  + " " + f(ncens.get(i)) + " " + f(nrisk.get(i)) + " "  + f(weight.get(i)) + " " + f(surv.get(i)) + " " + (varhaz.get(i)) + "  " + stderr.get(i) + "  " + stdlow.get(i) + "  " + lower.get(i) + "  " + upper.get(i) +"\r\n";
            
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
    public LinkedHashMap<Double,Integer> getNdead() {
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
