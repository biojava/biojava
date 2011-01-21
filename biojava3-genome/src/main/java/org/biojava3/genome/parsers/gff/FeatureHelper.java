/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.genome.parsers.gff;

import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FeatureHelper {

    /**
     * Build a list of individual features to allow easy indexing and to avoid iterating through large genome gff3 files
     * The index for the returned HashMap is the value of the attribute used to build the index
     * @param attribute
     * @param list
     * @return
     */
    static public LinkedHashMap<String,FeatureList> buildFeatureAtrributeIndex(String attribute,FeatureList list){

        LinkedHashMap<String,FeatureList> featureHashMap = new LinkedHashMap<String,FeatureList>();
        FeatureList featureList = list.selectByAttribute(attribute);
        for(FeatureI feature : featureList){
            String value = feature.getAttribute(attribute);
            FeatureList features = featureHashMap.get(value);
            if(features == null){
                features = new FeatureList();
                featureHashMap.put(value, features);
            }
            features.add(feature);
        }

        return featureHashMap;
    }

}
