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

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * Creates objects and returns them, and stores them in an internal
 * map of singletons for reference. Takes up a lot of memory!
 * @author Richard Holland
 * @since 1.5
 */
public class SimpleRichObjectBuilder implements RichObjectBuilder {
    
    private static Map objects = new HashMap();
    
    /**
     * {@inheritDoc}
     * Instantiates and returns objects, that's all there is to it.
     */
    public Object buildObject(Class clazz, List paramsList) {
        // put the class into the hashmap if not there already
        if (!objects.containsKey(clazz)) objects.put(clazz,new HashMap());
        Map contents = (Map)objects.get(clazz);
        // convert the params list to remove nulls as we can't process those.
        List ourParamsList = new ArrayList(paramsList);
        for (Iterator i = ourParamsList.iterator(); i.hasNext(); ) 
        	if (i.next()==null) i.remove();
        // return the constructed object from the hashmap if there already
        if (contents.containsKey(ourParamsList)) return contents.get(ourParamsList);
        // otherwise build it.
        try {
            // Load the class
            Class[] types = new Class[ourParamsList.size()];
            // Find its constructor with given params
            for (int i = 0; i < ourParamsList.size(); i++) {
                if (ourParamsList.get(i) instanceof Set) types[i] = Set.class;
                else if (ourParamsList.get(i) instanceof Map) types[i] = Map.class;
                else if (ourParamsList.get(i) instanceof List) types[i] = List.class;
                else types[i] = ourParamsList.get(i).getClass();
            }
            Constructor c = clazz.getConstructor(types);
            // Instantiate it with the parameters
            Object o = c.newInstance(ourParamsList.toArray());
            // store it for later in the singleton map
            contents.put(ourParamsList, o);
            // return it
            return o;
        } catch (Exception e) {
            StringBuffer paramsstuff = new StringBuffer();
            paramsstuff.append(clazz);
            paramsstuff.append("(");
            for (int i = 0; i < ourParamsList.size(); i++) {
                if (i>0) paramsstuff.append(",");
            	paramsstuff.append(ourParamsList.get(i).getClass());
            }
            paramsstuff.append(")");
            IllegalArgumentException ie = new IllegalArgumentException("Could not find constructor for "+paramsstuff);
            ie.initCause(e);
            throw ie;
        }
    }
    
}
