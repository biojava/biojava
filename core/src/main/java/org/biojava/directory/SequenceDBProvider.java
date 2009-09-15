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

import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.SequenceDBLite;

/**
 * Interfaces for named resources that can provide sequences via a
 * database given some configuration information as defined by the
 * OBDA standard.
 *
 * @author Thomas Down
 * @author Keith James
 * @author Matthew Pocock
 */
public interface SequenceDBProvider {
    /**
     * The name of this provider.
     *
     * @return the provider's name.
     */
    public String getName();
    
    /**
     * Get a sequence database.
     *
     * @param config a Map containing key-value pairs identifying the
     * database to resolve.
     * @return a SequenceDBLite that was resolved.
     */
    public SequenceDBLite getSequenceDB(Map config)
        throws RegistryException, BioException;
}
