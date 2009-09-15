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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.NoSuchElementException;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.program.indexdb.BioStore;
import org.biojava.bio.program.indexdb.Record;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.StreamReader;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;
import org.biojava.utils.io.RAF;
import org.biojava.utils.lsid.LifeScienceIdentifier;
import org.biojava.utils.lsid.LifeScienceIdentifierParseException;

/**
 * <code>FlatSequenceDB</code> is an OBDA flatfile sequence databank
 * implementation. It is backed by an index created using the
 * <code>org.biojava.bio.program.indexdb</code> package.
 *
 * @author Keith James
 */
public class FlatSequenceDB extends Unchangeable implements SequenceDBLite
{
    private BioStore index;
    private String dbName;
    private LifeScienceIdentifier format;

    public FlatSequenceDB(String location, String dbName)
        throws IOException, BioException
    {
        this.dbName = dbName;
        index = new BioStore(new File(location), false);

        try
        {
            Annotation config = index.getMetaData();
            String lsid = (String) config.getProperty("format");
            format = LifeScienceIdentifier.valueOf(lsid);
        }
        catch (NoSuchElementException nsee)
        {
            throw new BioException("Malformed OBDA index '"
                                   + location
                                   + "' does not indicate sequence format",nsee);
        }
        catch (LifeScienceIdentifierParseException lse)
        {
            throw new BioException("Malformed OBDA index '"
                                   + location
                                   + "' has a format identifier which is not a valid LSID",lse);
        }
    }

    public String getName()
    {
        return dbName;
    }

    public Sequence getSequence(String id)
        throws IllegalIDException, BioException
    {
        try
        {
            Record record = index.get(id);
            RAF seqRAF = record.getFile();
            int recLength = record.getLength();
            seqRAF.seek(record.getOffset());

            byte [] bytes = new byte [recLength];
            seqRAF.readFully(bytes, 0, recLength);
            InputStream is = new ByteArrayInputStream(bytes);

            int formatId = SeqIOTools.identifyFormat(format.getNamespaceId(),
                                                     format.getObjectId());

            SequenceFormat sf = SeqIOTools.getSequenceFormat(formatId);
            Alphabet alpha = SeqIOTools.getAlphabet(formatId);
            SymbolTokenization toke = alpha.getTokenization("token");
            SequenceBuilderFactory sbf = SeqIOTools.getBuilderFactory(formatId);

            SequenceIterator si = new StreamReader(is, sf, toke, sbf);
            return si.nextSequence();
        }
        catch (NoSuchElementException nsee)
        {
            throw new IllegalIDException("Failed to find sequence with ID "
                                         + id
                                         + " in database "
                                         + getName());
        }
        catch (IOException ioe)
        {
            throw new BioException("Failed to retrieve sequence with ID "
                                   + id, ioe);
        }
    }

    /**
     * <code>addSequence</code> always throws a
     * <code>ChangeVetoException</code> as this implementation is
     * immutable.
     *
     * @param sequence a <code>Sequence</code>.
     *
     * @exception ChangeVetoException
     */
    public void addSequence(Sequence sequence) throws ChangeVetoException
    {
        throw new ChangeVetoException("Failed to add sequence."
                                      + " Sequences may not be added"
                                      + " to a flat database");
    }

    /**
     * <code>removeSequence</code> always throws a
     * <code>ChangeVetoException</code> as this implementation is
     * immutable.
     *
     * @param id a <code>String</code>.
     *
     * @exception ChangeVetoException
     */
    public void removeSequence(String id) throws ChangeVetoException
    {
        throw new ChangeVetoException("Failed to add sequence."
                                      + " Sequences may not be removed"
                                      + " from a flat database");
    }
}
