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

package org.biojava.nbio.genome.parsers.gff;

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
