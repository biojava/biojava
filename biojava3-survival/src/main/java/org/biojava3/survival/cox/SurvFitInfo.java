/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox;


import java.util.LinkedHashMap;

/**
 * Contains info for graphing km figures
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SurvFitInfo {

    private LinkedHashMap<String, StrataInfo> strataInfoHashMap = new LinkedHashMap<String, StrataInfo>();

    /**
     * @return the strataInfoHashMap
     */
    public LinkedHashMap<String, StrataInfo> getStrataInfoHashMap() {
        return strataInfoHashMap;
    }

    /**
     * @param strataInfoHashMap the strataInfoHashMap to set
     */
    public void setStrataInfoHashMap(LinkedHashMap<String, StrataInfo> strataInfoHashMap) {
        this.strataInfoHashMap = strataInfoHashMap;
    }

    /**
     *
     * @param siHashMap
     * @param label
     */
    public void addStrataInfoHashMap(LinkedHashMap<String, StrataInfo> siHashMap,String label) {
        for (String key : siHashMap.keySet()) {
            StrataInfo si = siHashMap.get(key);
            strataInfoHashMap.put(label + " " + key, si);
        }
    }
}
