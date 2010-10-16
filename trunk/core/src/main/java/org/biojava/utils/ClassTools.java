/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of either the BSD licence or the GNU Lesser General
 * Public Licence.  These should be distributed with the code. 
 * If you do not have copies see:
 *
 *      http://www.opensource.org/licenses/bsd-license.php
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
 
package org.biojava.utils;

/**
 * Utility methods for manipulating class objects and resources.
 *
 * @author Thomas Down
 * @since 1.4
 */

public class ClassTools {
    private ClassTools() {
    }
    
    /**
     * Get the classloader which loaded the class of <code>obj</code>.
     */
    
    public static ClassLoader getClassLoader(Object obj) {
        return getClassLoader(obj.getClass());
    }
    
    /**
     * Get the classloader which loaded <code>clazz</code>.  This is
     * a "safe" method which handles <code>null</code> classloaders
     * and returns the system classloader instead.
     */
    
    public static ClassLoader getClassLoader(Class clazz) {
        ClassLoader cl = clazz.getClassLoader();
        if (cl == null) {
            cl = ClassLoader.getSystemClassLoader();
        }
        return cl;
    }
}
