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

package org.biojava.bio.program.ssbind;

import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.db.SequenceDBInstallation;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.bio.seq.impl.ViewSequence;

/**
 * <code>ViewSequenceFactory</code> is a base class for creating
 * search handlers which create and cache views on the query and
 * subject sequences.
 *
 * @author Keith James
 */
public abstract class ViewSequenceFactory
{
    // Supplier of instances of searched databases
    protected SequenceDBInstallation subjectDBs;
    // Holder for all query sequences
    protected SequenceDB querySeqHolder;
    // Identifier for database
    protected String databaseID;

    // Cache which holds view(s) of query sequence(s) which have
    // been instantiated for annotation
    protected Map queryViewCache;
    // Cache which holds view(s) of subject sequence(s) which have
    // been instantiated for annotation
    protected Map subjectViewCache;

    /**
     * <code>getQuerySeqHolder</code> returns the database of query
     * sequences used to retrieve sequences for creation of the
     * various result objects.
     *
     * @return a <code>SequenceDB</code> value.
     */
    public SequenceDB getQuerySeqHolder()
    {
        return querySeqHolder;
    }

    /**
     * <code>setQuerySeqHolder</code> sets the query sequence holder
     * to a specific database.
     *
     * @param querySeqHolder a <code>SequenceDB</code> containing the
     * query sequence(s).
     */
    public void setQuerySeqHolder(SequenceDB querySeqHolder)
    {
        this.querySeqHolder = querySeqHolder;
    }

    /**
     * <code>getSubjectDBInstallation</code> returns the installation
     * in which all the databases searched may be
     * found. <code>SequenceDB</code>s are retrieved for creation of
     * the various result objects.
     *
     * @return a <code>SequenceDBInstallation</code> containing the
     * subject database(s).
     */
    public SequenceDBInstallation getSubjectDBInstallation()
    {
        return subjectDBs;
    }

    /**
     * <code>setSubjectDBInstallation</code> sets the subject database
     * holder to a specific installation.
     *
     * @param subjectDBs a <code>SequenceDBInstallation</code>
     * containing the subject database(s)
     */
    public void setSubjectDBInstallation(SequenceDBInstallation subjectDBs)
    {
        this.subjectDBs = subjectDBs;
    }

    public void setDatabaseID(String databaseID)
    {
        this.databaseID = databaseID;
    }

    protected Sequence makeQueryViewSequence(String queryID)
        throws BioException
    {
        if (querySeqHolder == null)
            throw new BioException("Running with null query SequenceDB");

        if (queryViewCache.containsKey(queryID))
        {
            return (Sequence) queryViewCache.get(queryID);
        }
        else
        {
            Sequence query = querySeqHolder.getSequence(queryID);

            // It shouldn't happen, but it can with some implementations
            // of SequenceDB
            if (query == null)
                throw new BioException("Failed to retrieve query sequence "
                                       + "from holder using ID '"
                                       + queryID
                                       + "' (sequence was null)");

            query = new ViewSequence(query);
            queryViewCache.put(queryID, query);

            return query;
        }
    }

    protected Sequence makeSubjectViewSequence(String subjectID)
        throws BioException
    {
        if (subjectDBs == null)
            throw new BioException("Running with null subject "
                                   + "SequenceDB installation");

        SequenceDBLite subjectDB = subjectDBs.getSequenceDB(databaseID);

        // It shouldn't happen, but it can with some implementations
        // of SequenceDBInstallation
        if (subjectDB == null)
            throw new BioException("Failed to retrieve database "
                                   + "from installation using ID '"
                                   + databaseID
                                   + "' (database was null)");

        if (subjectViewCache.containsKey(subjectID))
        {
            return (Sequence) subjectViewCache.get(subjectID);
        }
        else
        {
            Sequence subject = subjectDB.getSequence(subjectID);

            // It shouldn't happen, but it can with some implementations
            // of SequenceDB
            if (subject == null)
                throw new BioException("Failed to retrieve subject sequence "
                                       + "from subjectDB using ID '"
                                       + subjectID
                                       + "' (sequence was null)");

            subject = new ViewSequence(subject);
            subjectViewCache.put(subjectID, subject);

            return subject;
        }
    }
}
