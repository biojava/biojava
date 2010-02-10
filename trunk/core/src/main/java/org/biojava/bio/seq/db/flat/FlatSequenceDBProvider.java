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

package org.biojava.bio.seq.db.flat;

import java.io.IOException;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.directory.RegistryException;
import org.biojava.directory.SequenceDBProvider;

/**
 * <code>FlatSequenceDBProvider</code> directory-services plugin for
 * flatfile databases.
 *
 * This class is instantiated automatically by the
 *                directory-services code, and is not of direct
 *                interest to users.
 *
 * @author Keith James
 */
public class FlatSequenceDBProvider implements SequenceDBProvider
{
    public String getName()
    {
        return "flat";
    }

    public SequenceDBLite getSequenceDB(Map config)
        throws RegistryException, BioException
    {
        String location = (String) config.get("location");
        if (location == null)
        {
            throw new RegistryException("Flat provider requires a"
                                        + " location parameter");
        }

        String dbName = (String) config.get("dbname");
        if (dbName == null)
        {
            throw new RegistryException("Flat provider requires a"
                                        + " dbname parameter");
        }

        try
        {
            return new FlatSequenceDB(location, dbName);
        }
        catch (IOException ioe)
        {
            throw new BioException("Flat provider failed to open index",ioe);
        }
    }
}
