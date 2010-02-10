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
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

/**
 * <p>A registry that loads up the standard biodirectory files.</p>
 *
 * <p>
 * This class will search for the following files in turn:
 * <ol>
 * <li>~/.bioinformatics/seqdatabase.ini where ~ is the JAVA user home system
 * property</li>
 * <li>/etc/bioinformatics/seqdatabase.ini</li>
 * <li>"http://www.open-bio.org/registry/seqdatabase.ini</li>
 * </ol>
 * </p>
 *
 * <p>The default search path may be replaced by an alternative search
 * path specified by the <code>OBDA_SEARCH_PATH</code> system
 * environment variable. This environment variable is a "+" delimted
 * string of files and URLs. The search order proceeds from read left
 * to right.</p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Keith James
 */
public class SystemRegistry {

    public static final String CONFIG_LOCATOR =
        "http://www.open-bio.org/registry/seqdatabase.ini";

    public static final String CONFIG_FILE = "seqdatabase.ini";

    public static final String OBDA_SEARCH_ENV = "OBDA_SEARCH_PATH";

    private static Registry systemRegistry;

    /**
     * Get the singleton Registry instance representing the system-wide
     * default registry.
     *
     * @return the system-wide Registry object.
     */
    public static Registry instance() {
        if (systemRegistry == null) {
            RegistryConfiguration.Composite regConfig
                = new RegistryConfiguration.Composite();
            Iterator i = getRegistryPath().iterator();

            while (i.hasNext()) {
                try {
                    String locator = (String) i.next();
                    URL url = new URL(locator);

                    if (url.getProtocol().equals("file")) {
                        File file = new File(url.getPath());
                        if (! file.exists() || ! file.canRead()) {
                            // FIXME - log this
                            continue;
                        }
                    }

                    BufferedReader stream = null;

                    try {
                        stream = new BufferedReader(new InputStreamReader(url.openStream()));
                    }
                    catch (IOException ioe) {
                        // FIXME - log this
                    }

                    if (stream != null) {
                        try {
                            RegistryConfiguration cfg =
                                OBDARegistryParser.parseRegistry(stream, locator);
                            regConfig.addBottomConfig(cfg);
                        } catch (Exception ex) {
                            // FIXME - log this
                            ex.printStackTrace();
                        }
                    } 
                } catch (Exception ex) {
                    // FIXME - log this
                    ex.printStackTrace();
                }
            }

            systemRegistry = new Registry(regConfig);
        }

        return systemRegistry;
    }

    /**
     * Get the list of places that will be searched for registry
     * files.
     *
     * @return a List of strings that are URLs to bioregistry files.
     */
    public static List getRegistryPath() {
        List registryPath = new ArrayList();

        String customPath = System.getProperty(OBDA_SEARCH_ENV);
        if (customPath != null) {
            StringTokenizer st = new StringTokenizer(customPath, "+");
            while (st.hasMoreTokens()) {
                registryPath.add(st.nextToken());
            }
        } else {
            String userHome = System.getProperty("user.home");
            if (userHome != null) {
                registryPath.add("file://"
                                 + userHome
                                 + "/.bioinformatics/"
                                 + CONFIG_FILE);
            }

            registryPath.add("file:///etc/bioinformatics/" + CONFIG_FILE);
            registryPath.add(CONFIG_LOCATOR);
        }

        return registryPath;
    }
}
