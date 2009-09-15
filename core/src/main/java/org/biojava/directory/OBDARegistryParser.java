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

package org.biojava.directory;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

/**
 * <p>This class encapsulates all the parsing of the OBDA registry
 * configuration file.</p>
 *
 * @author Brian Gilman
 * @author Thomas Down
 * @author Keith James
 * @version $Revision$
 */
public class OBDARegistryParser {

    /**
     * <code>parseRegistry</code> parses an Open Bioinformatics Database
     * Access (OBDA) configuration file.
     *
     * @param in a <code>BufferedReader</code>.
     * @param locator a <code>String</code> a configuration file
     * locator.
     *
     * @return a <code>RegistryConfiguration</code>.
     *
     * @exception IOException if an error reading the configuration
     * file.
     * @exception RegistryException if the configuration setup fails.
     */
    public static RegistryConfiguration parseRegistry(BufferedReader in,
                                                      String locator)
        throws IOException, RegistryException {
        String line = "";
        String dbName = "";
        String key = "";
        String value = "";
        Map config = new HashMap();
        Map currentDB = null;
	
        while ((line = in.readLine()) != null) {
            if (line.trim().length() > 0) {
                // We currently don't do anything with the version
                // number
                if (line.startsWith("VERSION=")) {
                    continue;
                }

                if (line.indexOf("[") > -1) {
                    dbName = line.substring(1, line.indexOf("]"));
                    currentDB = new HashMap();
                    currentDB.put("dbname", dbName);

                    // Create a List of Maps through which we can fall
                    // back if the first does not work
                    if (config.containsKey(dbName)) {
                        ((List) config.get(dbName)).add(currentDB);
                    } else {
                        List fallbacks = new ArrayList();
                        fallbacks.add(currentDB);
                        config.put(dbName, fallbacks);
                    }

                } else {
                    StringTokenizer strTok = new StringTokenizer(line, "=");
                    // Here we assume that there are only key = value
                    // pairs in the config file
                    key = strTok.nextToken();
                    if (strTok.hasMoreTokens()) {
                        value = strTok.nextToken();
                    } else {
                        value = "";
                    }

                    currentDB.put(key.trim(), value.trim());
                }
            }
        }

        return new RegistryConfiguration.Impl(locator, Collections.unmodifiableMap(config));
    }
}
