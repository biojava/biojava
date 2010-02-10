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

package org.biojava.bio.seq.db.biofetch;

import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.directory.RegistryException;
import org.biojava.directory.SequenceDBProvider;

/**
 * Directory-services plugin for biofetch databases.
 *
 * This class is instantiated automatically by the
 *                directory-services code, and is not of direct
 *                interest to users.
 *
 * @author Thomas Down
 * @author Keith James
 * @since 1.3
 */
public class BioFetchSequenceDBProvider implements SequenceDBProvider {
    public String getName() {
        return "biofetch";
    }

    public SequenceDBLite getSequenceDB(Map config)
        throws RegistryException, BioException
    {
        String location = (String) config.get("location");
        if (location == null) {
            throw new RegistryException("BioFetch provider requires"
                                        + " a location parameter");
        }

        String dbName = (String) config.get("dbname");
        if (dbName == null) {
            throw new RegistryException("BioFetch provider requires"
                                        + " a dbname parameter");
        }

        return new BioFetchSequenceDB(location, dbName);
    }
}
