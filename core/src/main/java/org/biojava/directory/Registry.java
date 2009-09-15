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

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.utils.ClassTools;
import org.biojava.utils.Services;

/**
 * <p><code>Registry</code> is a factory which gets implementations of
 * the BioJava <code>SequenceDBLite</code> interface. This is the
 * point of entry for OBDA access.</p>
 *
 * @author Brian Gilman
 * @author Thomas Down
 * @author Keith James
 *
 * @version $Revision$
 */
public class Registry {

    /**
     * Registry Configuration instance
     */
    private RegistryConfiguration regConfig = null;

    /**
     * Creates a new OBDA <code>Registry</code> with the specified
     * configuration.
     *
     * @param regConfig a <code>RegistryConfiguration</code>.
     */
    public Registry(RegistryConfiguration regConfig) {
        this.regConfig = regConfig;
    }

    /**
     * <code>getDatabase</code> retrieves a database instance known by
     * a name <code>String</code>.
     *
     * @param dbName a <code>String</code> database name.
     *
     * @return a <code>SequenceDBLite</code>.
     *
     * @exception RegistryException if the registry does not contain a
     * configuration for the specified name.
     * @exception BioException if the provider fails.
     */
    public SequenceDBLite getDatabase(String dbName)
        throws RegistryException, BioException {

        String providerName = "";

        List dbConfigs =
            (List) getRegistryConfiguration().getConfiguration().get(dbName);

        if (dbConfigs == null) {
            throw new RegistryException("Failed to find a configuration"
                                        + " for database: "
                                        + dbName);
        }

        for (Iterator ci = dbConfigs.iterator(); ci.hasNext();) {
            Map dbConfig = (Map) ci.next();
            providerName = (String) dbConfig.get("protocol");

            SequenceDBLite db = null;
            try {
                db = getProvider(providerName).getSequenceDB(dbConfig);
            } catch (RegistryException re) {
                // We allow RegistryExceptions to cause a fallback to
                // an alternative provider in the same config
                continue;
            }
            catch (Exception e) {
                // But more serious exceptions cause a failure
                throw new RegistryException("Failed to configure database "
                                            + dbName);
            }

            if (db != null)
                return db;
        }

        throw new RegistryException("Failed to find a configuration"
                                    + " for database: "
                                    + dbName);
    }

    private SequenceDBProvider getProvider(String providerName)
        throws RegistryException {
        try {
            ClassLoader loader = ClassTools.getClassLoader(this);
            Iterator implNames =
                Services.getImplementationNames(SequenceDBProvider.class, loader).iterator();

            while (implNames.hasNext()) {
              String className = (String) implNames.next();
              try {
                Class clazz = loader.loadClass(className);
                SequenceDBProvider seqDB =
                    (SequenceDBProvider) clazz.newInstance();
                if (seqDB.getName().equals(providerName)) {
                    return seqDB;
                }
              } catch (ClassNotFoundException ce) {
                throw new RegistryException(
                  "Could not find class: " + className +
                  " for service provider " + providerName, ce
                );
              }
            }

            throw new ProviderNotFoundException("No such provider exists: "
                                                + providerName);
        } catch (Exception e) {
            throw new RegistryException("Error accessing"
                                        + " SequenceDBProvider services",e);
        }
    }

    /**
     * <code>getRegistryConfiguration</code> returns the configuration
     * of the registry.
     *
     * @return a <code>RegistryConfiguration</code>.
     */
    public RegistryConfiguration getRegistryConfiguration() {
        return this.regConfig;
    }
}
