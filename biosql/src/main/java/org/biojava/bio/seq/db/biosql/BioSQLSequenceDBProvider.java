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

package org.biojava.bio.seq.db.biosql;

import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.directory.RegistryException;
import org.biojava.directory.SequenceDBProvider;

/**
 * @author Thomas Down
 * @deprecated Use hibernate and org.biojavax.bio.db.*
 */
public class BioSQLSequenceDBProvider implements SequenceDBProvider {
    public String getName() {
        return "biosql";
    }

    public SequenceDBLite getSequenceDB(Map config)
        throws RegistryException, BioException
    {
	String location = (String) config.get("location");
	if (location == null) {
	    throw new RegistryException("BioSQL provider requires"
                                    + " a 'location' parameter");
	}

	String dbName = (String) config.get("dbname");
	if (dbName == null) {
	    throw new RegistryException("BioSQL provider requires"
                                    + " a 'dbname' parameter");
	}

	String biodatabase = (String) config.get("biodbname");
	if (biodatabase == null) {
	    throw new RegistryException("BioSQL provider requires"
                                    + " a 'biodbname' parameter");
	}

	String userName = (String) config.get("user");
	if (userName == null) {
	    throw new RegistryException("BioSQL provider requires"
                                    + " a 'user' parameter");
	}

	String password = (String) config.get("passwd");
	if (password == null) {
	    throw new RegistryException("BioSQL provider requires"
                                    + " a 'passwd' parameter");
	}

	String driver = (String) config.get("driver");
	if (driver== null) {
	    throw new RegistryException("BioSQL provider requires"
                                    + " a 'driver' parameter");
	}

	return new BioSQLSequenceDB(driver,
																location,
                                userName,
                                password,
                                biodatabase,
                                false);
    }
}
