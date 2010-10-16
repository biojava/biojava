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

package org.biojava.bio.program.indexdb;

import java.util.List;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.utils.io.RAF;

/**
 * <code>IndexStore</code> is an interface for indexing flatfiles
 * according to the OBDA specification. It represents a map of Record instances
 * by a primary ID and any number of Records associated with an ID in some
 * seccondary namespace.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
public interface IndexStore {

    /**
     * <code>get</code> returns a record specified by a primary
     * identifier.
     *
     * @param id a <code>String</code> primary ID.
     *
     * @return a <code>Record</code>.
     *
     * @exception BioException if an error occurs or if there is no Record
     *            associated with the id
     */
    public Record get(String id) throws BioException;

    /**
     * <code>get</code> returns a list of <code>Record</code>s by
     * searching against the primary identifiers if the namespace
     * argument is equal to the primary namespace or otherwise by
     * searching the secondary namespaces. The list of Record instances retuned
     * may be empty, but is never null.
     *
     * @param id a <code>String</code> primary ID.
     * @param namespace a <code>String</code>.
     *
     * @return a <code>List</code> of <code>Record</code>s.
     *
     * @exception BioException if an error occurs.
     */
    public List get(String id, String namespace) throws BioException;

    /**
     * <code>getMetaData</code> returns a data structure which
     * represents an OBDA "config.dat" flatfile indexing configuration
     * file.
     *
     * @return an <code>Annotation</code>.
     */
    public Annotation getMetaData();

    /**
     * <code>writeRecord</code> creates and writes a new
     * <code>Record</code>
     *
     * @param file a <code>RAF</code> file.
     * @param offset a <code>long</code> byte offset.
     * @param length an <code>int</code> byte record length.
     * @param id a <code>String</code> primary ID.
     * @param secIDs a <code>Map</code> of primary ID to a
     * <code>List</code> of secondary IDs.
     */
    public void writeRecord(RAF file,
                            long offset,
                            int length,
                            String id,
                            Map secIDs);
}
