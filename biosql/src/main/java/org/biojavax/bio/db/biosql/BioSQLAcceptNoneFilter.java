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

package org.biojavax.bio.db.biosql;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Collections;
import java.util.Map;

import org.biojava.bio.seq.Feature;

/**
 * The class that accepts no features.
 * <p>
 * Use the FeatureFilter.none member.
 *
 * @author Matthew Pocock
 * @author Richard Holland
 * @since 1.5
 */
public class BioSQLAcceptNoneFilter implements BioSQLFeatureFilter {
    private Method isNull;
    
    protected BioSQLAcceptNoneFilter() {
        super();
        try {
            // Lazy load the Restrictions class.
            Class restrictions = Class.forName("org.hibernate.criterion.Restrictions");
            // Lookup the isNull method
            this.isNull = restrictions.getMethod("isNull", new Class[]{String.class});
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        } catch (NoSuchMethodException e) {
            throw new RuntimeException(e);
        }
    }
    
    public Object asCriterion() {
        try {
            // All feature have parents, so checking for null means none will ever be returned.
            return this.isNull.invoke(null,new Object[]{"parent"});
        } catch (InvocationTargetException e) {
            throw new RuntimeException(e);
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }
    }
    
    public Map criterionAliasMap() {
        return Collections.EMPTY_MAP;
    }
    
    public boolean equals(Object o) {
        return o instanceof BioSQLAcceptNoneFilter;
    }
    
    public boolean accept(Feature f) { return false; }
    
    public int hashCode() {
        return 1;
    }
    
    public String toString() {
        return "None";
    }
}

