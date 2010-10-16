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

package org.biojava.utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Set;

/**
 * Utility methods for handling META-INF/services files
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.3
 */

public class Services {
    /**
     * Return a Set of names of implementations of the
     * given service interface in the classloader from 
     * which BioJava was loaded.
     */
    public static Set getImplementationNames(Class serviceIF)
        throws IOException
    {
        return getImplementationNames(serviceIF,
                                      ClassTools.getClassLoader(Services.class));
    }

    /**
     * Return a List of names of implementations of the
     * given service interface available in a given
     * classloader.
     */
    public static Set getImplementationNames(Class serviceIF, ClassLoader loader)
        throws IOException
    {
        String serviceName = serviceIF.getName();
        Enumeration serviceFiles = loader.getResources("META-INF/services/"
                                                       + serviceName);
        Set names = new HashSet();

        while (serviceFiles.hasMoreElements()) {
            URL serviceFile = (URL) serviceFiles.nextElement();
            BufferedReader serviceReader =
                new BufferedReader(new InputStreamReader(serviceFile.openStream()));
            String implName;

            while ((implName = serviceReader.readLine()) != null) {
              if(implName.length() > 0) {
                names.add(implName);
              }
            }
        }

        return Collections.unmodifiableSet(names);
    }
}
