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
import java.util.Date;
import java.util.List;
import java.util.Map ;
import java.util.HashMap ;


/** a comparator to sort DasSources 
 * @author Andreas Prlic, Thomas Down
 */


public abstract class DasSourceComparator
    implements Comparator
{ 

    private final String name ;
    private static final Map COMPS_BY_NAME;

    private static final int TWODAYS = 1000 * 60 * 60 * 24 * 2;

    public DasSourceComparator(String str) {
	//System.out.println("new dasSourceComparator " + str);
	name = str ;
    }
   
    public static final Comparator BY_ID = new DasSourceComparator("id") {
        protected Comparable getField(DasSource ds) {
            return ds.getId();
        }
    };    

    public static final Comparator BY_NICKNAME = new DasSourceComparator("nickname") {
        protected Comparable getField(DasSource ds) {
            return ds.getNickname();
        }
    };    
    
    public static final Comparator BY_STATUS = new DasSourceComparator("status"){
		protected Comparable getField(DasSource ds) {
			
			Date now = new Date();
			
			if (ds.getLeaseDate().getTime() < (now.getTime() - TWODAYS))
				return new Integer(0);
			return new Integer(1);
		}
	};
    
    
    public static final Comparator BY_REGISTER_DATE = new DasSourceComparator("registerdate") {
        protected Comparable getField(DasSource ds) {
            return ds.getRegisterDate();
        }
    };
    public static final Comparator BY_LEASE_DATE = new DasSourceComparator("leasedate") {
        protected Comparable getField(DasSource ds) {
            return ds.getLeaseDate();
        }
    };
    public static final Comparator BY_URL = new DasSourceComparator("url") {
        protected Comparable getField(DasSource ds) {
            return ds.getUrl();
        }
    };
    public static final Comparator BY_ADMIN_EMAIL = new DasSourceComparator("adminemail") {
        protected Comparable getField(DasSource ds) {
            return ds.getAdminemail();
        }
    };
    public static final Comparator BY_DESCRIPTION = new DasSourceComparator("description") {
        protected Comparable getField(DasSource ds) {
            return ds.getDescription();
        }
    };
    public static final Comparator BY_CAPABILITIES = new DasSourceComparator("capabilities") {
        protected Comparable getField(DasSource ds) {
            List<String> caps = ds.getValidCapabilities();
            //System.out.println("return="+caps.length);
            
            return caps==null ? 0 :caps.size();
        }
    };
    public static final Comparator BY_COORDINATE_SYSTEM = new DasSourceComparator("coordinateSystem") {
        protected Comparable getField(DasSource ds) {
            List<DasCoordinateSystem> dcss = ds.getCoordinateSystem();
            return dcss.size() == 0 ? "" : dcss.get(0).toString();
        }
    };

    static {
        COMPS_BY_NAME = new HashMap();
        COMPS_BY_NAME.put(BY_ID.toString(),                BY_ID);
        COMPS_BY_NAME.put(BY_NICKNAME.toString(),          BY_NICKNAME);
        COMPS_BY_NAME.put(BY_REGISTER_DATE.toString(),     BY_REGISTER_DATE);
        COMPS_BY_NAME.put(BY_LEASE_DATE.toString(),        BY_LEASE_DATE);
        COMPS_BY_NAME.put(BY_URL.toString(),               BY_URL);
        COMPS_BY_NAME.put(BY_ADMIN_EMAIL.toString(),       BY_ADMIN_EMAIL);
        COMPS_BY_NAME.put(BY_DESCRIPTION.toString(),       BY_DESCRIPTION);
        COMPS_BY_NAME.put(BY_CAPABILITIES.toString(),      BY_CAPABILITIES);
        COMPS_BY_NAME.put(BY_COORDINATE_SYSTEM.toString(), BY_COORDINATE_SYSTEM);
        COMPS_BY_NAME.put(BY_STATUS.toString(),            BY_STATUS);
    }

   

    public static Comparator fromString(String name) {
        if (COMPS_BY_NAME.containsKey(name)) {
            return (Comparator) COMPS_BY_NAME.get(name);
        } else {
            throw new IllegalArgumentException("Can't compare by key " + name);
        }
    }

    protected abstract Comparable getField(DasSource ds);

    /** compare two DasSource objects */
    public int compare( Object a, Object b) {
        
        DasSource x = (DasSource) a ;
        DasSource y = (DasSource) b ;
        return getField(x).compareTo(getField(y));
    }

    public String toString() {
        return name;
    }


}


