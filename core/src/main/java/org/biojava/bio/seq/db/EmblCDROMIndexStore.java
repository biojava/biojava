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

package org.biojava.bio.seq.db;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.emblcd.DivisionLkpReader;
import org.biojava.bio.seq.db.emblcd.EmblCDROMIndexReader;
import org.biojava.bio.seq.db.emblcd.EmblCDROMRandomAccess;
import org.biojava.bio.seq.db.emblcd.EntryNamIdxReader;
import org.biojava.bio.seq.db.emblcd.EntryNamRandomAccess;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.SymbolTokenization;

/**
 * <p><code>EmblCDROMIndexStore</code>s implement a read-only
 * <code>IndexStore</code> backed by EMBL CD-ROM format binary
 * indices. The required index files are typically named
 * "division.lkp" and "entrynam.idx". As an <code>IndexStore</code>
 * performs lookups by sequence ID, the index files "acnum.trg" and
 * "acnum.hit" (which store additional accession number data) are not
 * used.</p>
 *
 * <p>The sequence IDs are found using a binary search via a pointer
 * into the index file. The whole file is not read unless a request
 * for all the IDs is made using the getIDs() method. The set of IDs
 * is then cached after the first pass. This class also has a
 * <code>close()</code> method to free resources associated with the
 * underlying <code>RandomAccessFile</code>.</p>
 *
 * <p>The binary index files may be created using the EMBOSS programs
 * dbifasta, dbiblast, dbiflat or dbigcg. The least useful from the
 * BioJava perspective is dbigcg because we do not have a
 * <code>SequenceFormat</code> implementation for GCG format
 * files.</p>
 *
 * <p>The <code>Index</code> instances returned by this class do not
 * have the record length set because this information is not
 * available in the binary index. The value -1 is used instead, as
 * described in the <code>Index</code> interface.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public class EmblCDROMIndexStore implements IndexStore 
{
    private File divisionLkp;
    private File entryNamIdx;

    // Optional PATH prefix to append to the filename(s) extracted
    // from the binary indices
    private File pathPrefix;

    private SequenceFormat         format;
    private SequenceBuilderFactory factory;
    private SymbolTokenization     parser;

    // Maps the file numbers used in the indices to the real file names
    private Map seqFiles;
    // Set view of file names
    private Set fileSet;
    // Lazily instantiated if someone asks for all the IDs at once
    private Set seqIds;
    // The database name defined in the index header
    private String name;

    // Details of the master index records
    private long divRecordCount;
    // Details of the ID/offset records
    private int entryRecordLength;
    private long entryRecordCount;

    // The random access file containing ID/offset records
    private EmblCDROMRandomAccess entryRandomAccess;

    /**
     * Creates a new <code>EmblCDROMIndexStore</code> backed by a
     * random access binary index.
     *
     * @param divisionLkp a <code>File</code> containing the master
     * index.
     * @param entryNamIdx a <code>File</code> containing the sequence
     * IDs and offsets.
     * @param format a <code>SequenceFormat</code>.
     * @param factory a <code>SequenceBuilderFactory</code>.
     * @param parser a <code>SymbolTokenization</code>.
     *
     * @exception IOException if an error occurs.
     */
    public EmblCDROMIndexStore(File                   divisionLkp,
                               File                   entryNamIdx,
                               SequenceFormat         format,
                               SequenceBuilderFactory factory,
                               SymbolTokenization     parser)
        throws IOException
    {
        // Set to the empty abstract path
        this(new File(""), divisionLkp, entryNamIdx,
             format, factory, parser);
    }

    /**
     * Creates a new <code>EmblCDROMIndexStore</code> backed by a
     * random access binary index.
     *
     * @param pathPrefix a <code>File</code> containing the abstract
     * path to be appended to sequence database filenames retrieved
     * from the binary index.
     * @param divisionLkp a <code>File</code> containing the master
     * index.
     * @param entryNamIdx a <code>File</code> containing the sequence
     * IDs and offsets.
     * @param format a <code>SequenceFormat</code>.
     * @param factory a <code>SequenceBuilderFactory</code>.
     * @param parser a <code>SymbolTokenization</code>.
     *
     * @exception IOException if an error occurs.
     */
    public EmblCDROMIndexStore(File                   pathPrefix,
                               File                   divisionLkp,
                               File                   entryNamIdx,
                               SequenceFormat         format,
                               SequenceBuilderFactory factory,
                               SymbolTokenization     parser)
        throws IOException
    {
        this.divisionLkp = divisionLkp;
        this.entryNamIdx = entryNamIdx;
        this.format      = format;
        this.factory     = factory;
        this.parser      = parser;
        this.pathPrefix  = pathPrefix;

        initialise();
    }

    /**
     * <code>getPathPrefix</code> returns the abstract path currently
     * being appended to the raw sequence database filenames extracted
     * from the binary index. This value defaults to the empty
     * abstract path.
     *
     * @return a <code>File</code>.
     */
    public File getPathPrefix()
    {
        return pathPrefix;
    }

    /**
     * <code>setPathPrefix</code> sets the abstract path to be
     * appended to sequence database filenames retrieved from the
     * binary index. E.g. if the binary index refers to the database
     * as 'SWALL' and the <code>pathPrefix</code> is set to
     * "/usr/local/share/data/seq/", then the <code>IndexStore</code>
     * will know the database path as
     * "/usr/local/share/data/seq/swall" and any <code>Index</code>
     * instances produced by the store will return the latter path
     * when their getFile() method is called. This value defaults to
     * the empty abstract path.
     *
     * @param pathPrefix a <code>File</code> prefix specifying the
     * abstract path to append.
     */
    public void setPathPrefix(File pathPrefix)
    {
        this.pathPrefix = pathPrefix;
    }

    /**
     * <code>getName</code> returns the database name as defined
     * within the EMBL CD-ROM index.
     *
     * @return a <code>String</code> value.
     */
    public String getName()
    {
        return name;
    }

    /**
     * <code>store</code> adds an <code>Index</code> to the store. As
     * EMBL CD-ROM indices are read-only, this implementation throws a
     * <code>BioException</code>.
     *
     * @param index an <code>Index</code>.
     *
     * @exception IllegalIDException if an error occurs.
     * @exception BioException if an error occurs.
     */
    public void store(Index index)
        throws IllegalIDException, BioException
    {
        throw new BioException("Failed to add Index: store is read-only."
                               + " To add sequences use the dbi programs"
                               + " supplied in EMBOSS");
    }

    /**
     * <code>commit</code> commits changes. As EMBL CD-ROM indices are
     * read-only, this implementation throws a
     * <code>BioException</code>.
     *
     * @exception BioException if an error occurs.
     */
    public void commit() throws BioException
    {
        throw new BioException("Failed to commit: store is read-only."
                               + " To add sequences use the dbi programs"
                               + " supplied in EMBOSS");
    }

    /**
     * <code>rollback</code> rolls back changes made since the last
     * <code>commit</code>. As EMBL CD-ROM indices are read-only, this
     * implementation does nothing.
     */
    public void rollback() { }

    public Index fetch(String id) throws IllegalIDException, BioException
    {
        Index index = null;

        try
        {
            Object [] enRecord = entryRandomAccess.findRecord(id);

            if (enRecord.length == 0)
                throw new IllegalIDException("Failed to find ID: " + id);

            // Append current pathPrefix
            index =
                new SimpleIndex(new File(pathPrefix,
                                         (String) seqFiles.get((Integer)
                                                               enRecord[3])),
                                ((Long) enRecord[1]).longValue(), -1, id);
        }
        catch (IOException ioe)
        {
            throw new BioException("Failed to retrieve index for ID: " + id);
        }

        return index;
    }

    public Set getIDs()
    {
        if (seqIds == null)
        {
            seqIds = new HashSet((int) entryRecordCount);

            BufferedInputStream bis = null;

            try
            {
                bis =
                    new BufferedInputStream(new FileInputStream(entryNamIdx));
                EmblCDROMIndexReader ent = new EntryNamIdxReader(bis);

                for (long i = 0; i < entryRecordCount; i++)
                {
                    Object [] enRecord = ent.readRecord();
                    seqIds.add((String) enRecord[0]);
                }

                bis.close();
            }
            // File was not found, so don't try to close it
            catch (FileNotFoundException fnfe)
            {
                System.err.println("Failed to find file "
                                   + entryNamIdx.getName());
                fnfe.printStackTrace();
            }
            // File was opened, so try to close it
            catch (IOException ioe)
            {
                try
                {
                    bis.close();
                }
                catch (IOException ioe2)
                {
                    System.err.println("Failed to close input stream from file "
                                       + entryNamIdx.getName());
                }

                System.err.println("Failed to read file "
                                   + entryNamIdx.getName());
                ioe.printStackTrace();
            }
        }

        return Collections.unmodifiableSet(seqIds);
    }

    public Set getFiles()
    {
        return Collections.unmodifiableSet(fileSet);
    }

    public SequenceFormat getFormat()
    {
        return format;
    }

    public SequenceBuilderFactory getSBFactory()
    {
        return factory;
    }

    public SymbolTokenization getSymbolParser()
    {
        return parser;
    }

    /**
     * <code>close</code> closes the underlying
     * <code>EntryNamRandomAccess</code> which in turn closes the
     * lower level <code>RandomAccessFile</code>. This frees the
     * resources associated with the file.
     *
     * @exception IOException if an error occurs.
     */
    public void close() throws IOException
    {
        entryRandomAccess.close();
    }

    /**
     * <code>initialise</code> reads the headers of the index files to
     * obtain data about the record sizes and counts, database name
     * and sequence filenames. It then opens a random access file to
     * the ID index for lookups.
     *
     * @exception IOException if an error occurs.
     */
    private void initialise() throws IOException
    {
        BufferedInputStream bis = null;

        // First try to get details of file names and numbers from
        // master index file.
        try
        {
            bis = new BufferedInputStream(new FileInputStream(divisionLkp));
            EmblCDROMIndexReader div = new DivisionLkpReader(bis);

            divRecordCount  = div.readRecordCount();

            // The database name is the same in all the index headers
            name = div.readDBName();

            seqFiles = new HashMap((int) divRecordCount);

            // Store the file number->name mapping
            for (long i = divRecordCount; --i >= 0;)
            {
                Object [] divRecord = div.readRecord();

                Integer fileNumber = (Integer) divRecord[0];
                String    fileName = (String)  divRecord[1];

                seqFiles.put(fileNumber, fileName);
            }

            // Keep a Set view
            fileSet = new HashSet((int) divRecordCount);
            fileSet.addAll(seqFiles.values());

            bis.close();
        }
        // File was not found, so don't try to close it
        catch (FileNotFoundException fnfe)
        {
            System.err.println("Failed to find file "
                               + divisionLkp.getName());
            // Rethrow
            throw fnfe;
        }
        // File was opened, so try to close it
        catch (IOException ioe)
        {
            try
            {
                bis.close();
            }
            catch (IOException ioe2)
            {
                System.err.println("Failed to close input stream from file "
                                   + divisionLkp.getName());
            }

            System.err.println("Failed to read full set of sequence IDs file "
                               + divisionLkp.getName());
            // Rethrow
            throw ioe;
        }

        // Now try to get details of sequence ID index file
        try
        {
            bis = new BufferedInputStream(new FileInputStream(entryNamIdx));
            EmblCDROMIndexReader ent = new EntryNamIdxReader(bis);

            entryRecordLength = ent.readRecordLength();
            entryRecordCount  = ent.readRecordCount();

            bis.close();
        }
        // File was not found, so don't try to close it
        catch (FileNotFoundException fnfe)
        {
            System.err.println("Failed to find file "
                               + entryNamIdx.getName());
            // Rethrow
            throw fnfe;
        }
        // File was opened, so try to close it
        catch (IOException ioe)
        {
            try
            {
                bis.close();
            }
            catch (IOException ioe2)
            {
                System.err.println("Failed to close input stream from file "
                                   + entryNamIdx.getName());
            }

            System.err.println("Failed to read file "
                               + entryNamIdx.getName());
            // Rethrow
            throw ioe;
        }

        // Try to set up random access file
        try
        {
            entryRandomAccess = new EntryNamRandomAccess(entryNamIdx,
                                                         300,
                                                         entryRecordLength,
                                                         entryRecordCount);
        }
        // File was not found, so don't try to close it
        catch (FileNotFoundException fnfe)
        {
            System.err.println("Failed to find file "
                               + entryNamIdx.getName());
            try
            {
                bis.close();
            }
            catch (IOException ioe2)
            {
                System.err.println("Failed to close random access file "
                                   + entryNamIdx.getName());
            }
            // Rethrow
            throw fnfe;
        }
    }
}
