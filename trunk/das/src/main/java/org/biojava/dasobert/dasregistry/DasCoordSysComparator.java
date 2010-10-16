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
 * Created on 15.04.2004
 * @author Andreas Prlic
 *
 */


package org.biojava.dasobert.dasregistry ;

import java.util.Comparator ;
import java.util.Map ;
import java.util.HashMap ;

import org.biojava.dasobert.dasregistry.DasCoordinateSystem;

/** a comparator to sort DasSources 
 * @author Andreas Prlic 
 */


public abstract class DasCoordSysComparator
implements Comparator
{ 

	private final String name ;
	private static final Map COMPS_BY_NAME;


	public DasCoordSysComparator(String str) {
		//System.out.println("new dasSourceComparator " + str);
		name = str ;
	}

	
	
	
	public static final Comparator BY_NAME = new DasCoordSysComparator("name") {
		protected Comparable getField(DasCoordinateSystem ds) {
			return ds.getName();
		}
	};    

	public static final Comparator BY_ID = new DasCoordSysComparator("id") {
		protected Comparable getField(DasCoordinateSystem ds) {
			return ds.getUniqueId();
		}
	};    
	public static final Comparator BY_CATEGORY = new DasCoordSysComparator("category") {
		protected Comparable getField(DasCoordinateSystem ds) {
			return ds.getCategory();
		}
	};
	public static final Comparator BY_ORGANISM = new DasCoordSysComparator("organism") {
		protected Comparable getField(DasCoordinateSystem ds) {
			return ds.getOrganismName();
		}
	};
	public static final Comparator BY_TAXID = new DasCoordSysComparator("taxid") {
		protected Comparable getField(DasCoordinateSystem ds) {
			return ds.getNCBITaxId()+"";
		}
	};



	static {
		COMPS_BY_NAME = new HashMap();
		COMPS_BY_NAME.put(BY_ID.toString(),           BY_ID);
		COMPS_BY_NAME.put(BY_NAME.toString(),         BY_NAME);
		COMPS_BY_NAME.put(BY_CATEGORY.toString(),     BY_CATEGORY);
		COMPS_BY_NAME.put(BY_ORGANISM.toString(),     BY_ORGANISM);
		COMPS_BY_NAME.put(BY_TAXID.toString(),        BY_TAXID);
	}



	public static Comparator fromString(String name) {
		if (COMPS_BY_NAME.containsKey(name)) {
			return (Comparator) COMPS_BY_NAME.get(name);
		} else {
			throw new IllegalArgumentException("Can't compare by key " + name);
		}
	}

	protected abstract Comparable getField(DasCoordinateSystem ds);

	/** compare two DasCoordSys objects */
	public int compare( Object a, Object b) {
		DasCoordinateSystem x = (DasCoordinateSystem) a ;
		DasCoordinateSystem y = (DasCoordinateSystem) b ;
		return getField(x).compareTo(getField(y));
	}

	public String toString() {
		return name;
	}


}


