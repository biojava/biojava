/*
 *                  BioJava development code
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
 * Created on Dec 8, 2006
 * 
 */

package org.biojava.dasobert.dasregistry;

import java.util.ArrayList;
import java.util.List;

/** a class to filter a set of DAS sources, or to check if
 *  single DAS sources fulfill certain requirements
 * 
 * 
 * @author Andreas Prlic
 *
 */
public class DasSourceFilter {

	public boolean hasAuthority(DasSource source, String authority){
		if ( authority == null )
			return true;
		if ( authority.equals(""))
			return true;

		DasCoordinateSystem[] coords = source.getCoordinateSystem();
		for ( int j =0 ; j < coords.length ; j++ ){
			DasCoordinateSystem cs = coords[j];
			if ( authority.equalsIgnoreCase(cs.getName())) {
				return true;
			}
		}

		return false;
	}

	public boolean hasCapability(DasSource source, String capability){

		if ( capability == null)
			return true;
		if ( capability.equals(""))
			return true;

		return source.hasCapability(capability);
	}


	public boolean hasLabel(DasSource source, String label){
		if ( label == null)
			return true;
		if ( label.equals(""))
			return true;


		String[] labels = source.getLabels();
		for ( int j = 0 ; j< labels.length ; j++){
			String l = labels[j];
			if ( l.equalsIgnoreCase(label))
				return true;
		}

		return false;
	}
	
	public boolean hasType(DasSource source, String type){
		if ( type == null )
			return true;
		if ( type.equals(""))
			return true;

		DasCoordinateSystem[] coords = source.getCoordinateSystem();
		for ( int j =0 ; j < coords.length ; j++ ){
			DasCoordinateSystem cs = coords[j];
			if ( type.equalsIgnoreCase(cs.getCategory())) {
				return true;
			}
		}

		return false;
	}

	public boolean hasOrganism(DasSource source, String organism){
//		test for correct organism
		if ( organism == null)
			return true;
		if ( organism.equals(""))
			return true;


		DasCoordinateSystem[] coords = source.getCoordinateSystem();
		for ( int j =0 ; j < coords.length ; j++ ){
			DasCoordinateSystem cs = coords[j];
			if ( ( organism.equalsIgnoreCase(cs.getOrganismName())) ||
					( organism.equalsIgnoreCase(cs.getNCBITaxId()+""))) {
				return true;
			}
		}
		return false;
	}



	/** filter a set of DasSources by particular requirements
	 * all arguments can be set to null which means they are ignored
	 * 
	 * @param sources
	 * @param label
	 * @param organism
	 * @param authority
	 * @param capability
	 * @param type 
	 * @return an array of DasSources that match the requested filtering rules
	 */
	public  DasSource[] filterBy(DasSource[] sources, 
			String label, 
			String organism,
			String authority,
			String capability,
			String type) {

		if ( (label == null ) &&
				( organism == null) &&
				( authority == null) &&
				( capability == null) &&
				( type == null))
			return sources;


		List lst = new ArrayList();
		for (int i = 0 ; i < sources.length; i++) {
			DasSource source = sources[i];

			// test for correct label
			if (  hasLabel(source, label) &&
					hasOrganism(source, organism) &&
					hasAuthority(source,authority) &&
					hasCapability(source,capability) &&
					hasType(source,type)){
				lst.add(source);
			}
		}

		return (DasSource[]) lst.toArray(new DasSource[lst.size()]);

	}



}
