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

package org.biojavax;

import java.util.List;

/**
 * This interface allows a class to generate Rich objects based on a class
 * name and some parameters. The idea is to allow maintenance of a singleton
 * map.
 * @author Richard Holland
 * @since 1.5
 */
public interface RichObjectBuilder {
    
    /**
     * This method takes a class name and some parameters, and uses that
     * information to construct and return an equivalent object, usually by
     * calling the constructor on the class with the supplied parameters.
     * Note that it only works with classes whose constructors take only
     * Objects, and not primitives. It should return singletons.
     * @param clazz the class to instantiate and build
     * @param paramsList the parameters to pass to the constructor
     * @return an instance of the requested class/params combination. May or
     * may not be a singleton, but usually will be given as that is the 
     * purpose of this class.
     */
    public Object buildObject(Class clazz, List paramsList);
    
}
