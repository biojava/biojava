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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.program.tagvalue.ChangeTable;
import org.biojava.bio.program.tagvalue.Indexer;
import org.biojava.bio.program.tagvalue.LineSplitParser;
import org.biojava.bio.program.tagvalue.Parser;
import org.biojava.bio.program.tagvalue.ValueChanger;
import org.biojava.bio.seq.io.SeqIOConstants;
import org.biojava.utils.CommitFailure;
import org.biojava.utils.ParserException;
import org.biojava.utils.io.CountedBufferedReader;
import org.biojava.utils.io.RAF;
import org.biojava.utils.lsid.LifeScienceIdentifier;

/**
 * <code>IndexTools</code> contains static utility methods for
 * creating flatfile indices according to the OBDA standard.
 *
 * @author Keith James
 * @author Matthew Pocock
 */
public class IndexTools
{
    // Cannot be instantiated
    private IndexTools() { }

    /**
     * <code>indexFasta</code> indexes DNA, RNA or protein Fasta
     * format sequence files on primary identifier.
     *
     * @param location a <code>File</code> directory which will
     * contain the indices.
     * @param seqFiles a <code>File []</code> array of files to index.
     * @param alphabetIdentifier an <code>int</code> indicating the
     * type of sequence to be indexed. May be one of
     * <code>SeqIOConstants.DNA SeqIOConstants.RNA
     * SeqIOConstants.AA</code>.
     * @param name a <code>String</code> arbitrary database name.
     *
     * @exception FileNotFoundException if an error occurs.
     * @exception IOException if an error occurs.
     * @exception ParserException if an error occurs.
     * @exception BioException if an error occurs.
     */
    public static void indexFasta(String name, File location, File [] seqFiles,
                                  int alphabetIdentifier)
        throws FileNotFoundException, IOException, ParserException,
               BioException
    {
        BioStoreFactory bsf = new BioStoreFactory();
        bsf.setStoreName(name);

        switch (alphabetIdentifier)
        {
            case (SeqIOConstants.DNA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_FASTA_DNA);
                break;
            case (SeqIOConstants.RNA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_FASTA_RNA);
                break;
            case (SeqIOConstants.AA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_FASTA_AA);
                break;

            default:
                throw new IllegalArgumentException("Unknown alphabet identifier '"
                                                   + alphabetIdentifier
                                                   + "'");
        }

        _indexFasta(bsf, location, seqFiles);
    }

    /**
     * <code>indexEmbl</code> indexes DNA, RNA or protein EMBL format
     * sequence files on ID as primary identifier and AC as secondary.
     *
     * @param location a <code>File</code> directory which will
     * contain the indices.
     * @param seqFiles a <code>File []</code> array of files to index.
     * @param alphabetIdentifier an <code>int</code> indicating the
     * type of sequence to be indexed. May be one of
     * <code>SeqIOConstants.DNA SeqIOConstants.RNA
     * SeqIOConstants.AA</code>.
     * @param name a <code>String</code> arbitrary database name.
     *
     * @exception FileNotFoundException if an error occurs.
     * @exception IOException if an error occurs.
     * @exception ParserException if an error occurs.
     * @exception BioException if an error occurs.
     */
    public static void indexEmbl(String name, File location, File [] seqFiles,
                                 int alphabetIdentifier)
        throws FileNotFoundException, IOException, ParserException,
               BioException
    {
        BioStoreFactory bsf = new BioStoreFactory();
        bsf.setStoreName(name);

        switch (alphabetIdentifier)
        {
            case (SeqIOConstants.DNA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_EMBL_DNA);
                break;
            case (SeqIOConstants.RNA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_EMBL_RNA);
                break;
            case (SeqIOConstants.AA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_EMBL_AA);
                break;

            default:
                throw new IllegalArgumentException("Unknown alphabet identifier '"
                                                   + alphabetIdentifier
                                                   + "'");
        }

        _indexEmblLike(bsf, location, seqFiles);
    }

    /**
     * <code>indexGenbank</code> indexes DNA, RNA or protein Genbank
     * format sequence files on LOCUS as primary identifier and
     * ACCESSION as secondary.
     *
     * @param location a <code>File</code> directory which will
     * contain the indices.
     * @param seqFiles a <code>File []</code> array of files to index.
     * @param alphabetIdentifier an <code>int</code> indicating the
     * type of sequence to be indexed. May be one of
     * <code>SeqIOConstants.DNA SeqIOConstants.RNA
     * SeqIOConstants.AA</code>.
     * @param name a <code>String</code> arbitrary database name.
     *
     * @exception FileNotFoundException if an error occurs.
     * @exception IOException if an error occurs.
     * @exception ParserException if an error occurs.
     * @exception BioException if an error occurs.
     */
    public static void indexGenbank(String name, File location, File [] seqFiles,
                                    int alphabetIdentifier)
        throws FileNotFoundException, IOException, ParserException,
               BioException
    {
        BioStoreFactory bsf = new BioStoreFactory();
        bsf.setStoreName(name);

        switch (alphabetIdentifier)
        {
            case (SeqIOConstants.DNA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_GENBANK_DNA);
                break;
            case (SeqIOConstants.RNA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_GENBANK_RNA);
                break;
            case (SeqIOConstants.AA):
                bsf.setSequenceFormat(SeqIOConstants.LSID_GENBANK_AA);
                break;

            default:
                throw new IllegalArgumentException("Unknown alphabet identifier '"
                                                   + alphabetIdentifier
                                                   + "'");
        }

        _indexGenbank(bsf, location, seqFiles);
    }


    /**
     * <code>indexSwissprot</code> indexes Swissprot format protein
     * sequence files on ID as primary identifier.
     *
     * @param location a <code>File</code> directory which will
     * contain the indices.
     * @param seqFiles a <code>File []</code> array of files to index.
     * @exception FileNotFoundException if an error occurs.
     * @exception IOException if an error occurs.
     * @exception ParserException if an error occurs.
     * @exception BioException if an error occurs.
     */
    public static void indexSwissprot(String name, File location, File [] seqFiles)
        throws FileNotFoundException, IOException, ParserException,
               BioException
    {
        BioStoreFactory bsf = new BioStoreFactory();
        bsf.setStoreName(name);
        bsf.setSequenceFormat(LifeScienceIdentifier.valueOf("open-bio.org",
                                                            "swiss",
                                                            "protein" ));
        _indexEmblLike(bsf, location, seqFiles);
    }

    private static void _indexFasta(BioStoreFactory bsf,
                                    File location, File [] seqFiles)
       throws FileNotFoundException, IOException, BioException
    {
        bsf.setPrimaryKey("ID");
        bsf.setStoreLocation(location);
        bsf.addKey("ID", 10);

        BioStore store = bsf.createBioStore();

        for (int i = 0; i < seqFiles.length; i++)
        {
            // File data
            long newOffset = 0L;
            long oldOffset = 0L;
            RAF raf = new RAF(seqFiles[i], "r");
            Map map = new HashMap();

            CountedBufferedReader reader =
                new CountedBufferedReader(new FileReader(raf.getFile()));

            // Record data
            String id = "";

            String line = null;
            while ((line = reader.readLine()) != null)
            {
                if (line.startsWith(">"))
                {
                    // Write at end of record
                    if (newOffset > 0)
                    {
                        store.writeRecord(raf, oldOffset,
                                          (int) (newOffset - oldOffset),
                                          id, map);
                        oldOffset = newOffset;
                    }
                    newOffset = reader.getFilePointer();

                    int delimeter = line.indexOf(" ");
                    if (delimeter < 0)
                        id = line.substring(1);
                    else
                        id = line.substring(1, delimeter);
                }
                else
                {
                    newOffset = reader.getFilePointer();
                }
            }

            // Write final record
            store.writeRecord(raf, oldOffset,
                              (int) (newOffset - oldOffset),
                              id, map);
        }

        try
        {
            store.commit();
        }
        catch (CommitFailure ne)
        {
            throw new BioException("Failed to commit new index to file",ne);
        }
    }

    private static void _indexEmblLike(BioStoreFactory bsf,
                                       File location, File [] seqFiles)
        throws FileNotFoundException, IOException, ParserException,
               BioException
    {
        bsf.setPrimaryKey("ID");
        bsf.setStoreLocation(location);
        bsf.addKey("AC", 10);
        bsf.addKey("ID", 10);

        BioStore store = bsf.createBioStore();

        for (int i = 0; i < seqFiles.length; i++)
        {
            Indexer indexer = new Indexer(seqFiles[i], store);
            indexer.setPrimaryKeyName("ID");
            indexer.addSecondaryKey("AC");

            ChangeTable changeTable = new ChangeTable();

            changeTable.setChanger("ID", new ChangeTable.Changer()
                {
                    public Object change(Object value)
                    {
                        String s = (String) value;
                        int i = s.indexOf(" ");

                        if (i < 0)
                            return s;
                        else
                            return s.substring(0, i);
                    }
                });

            changeTable.setChanger("AC", new ChangeTable.Changer()
                {
                    public Object change(Object value)
                    {
                        String s = (String) value;
                        int i = s.indexOf(";");
                        return s.substring(0, i);
                    }
                });

            ValueChanger changer = new ValueChanger(indexer, changeTable);
            Parser parser = new Parser();

            while(parser.read(indexer.getReader(),
                              LineSplitParser.EMBL, changer));
        }

        try
        {
            store.commit();
        }
        catch (CommitFailure ne)
        {
            throw new BioException("Failed to commit new index to file",ne);
        }
    }

    private static void _indexGenbank(BioStoreFactory bsf,
                                      File location, File [] seqFiles)
        throws FileNotFoundException, IOException, ParserException,
               BioException
    {
        bsf.setPrimaryKey("LOCUS");
        bsf.setStoreLocation(location);
        bsf.addKey("LOCUS", 10);
        bsf.addKey("ACCESSION", 10);

        BioStore store = bsf.createBioStore();

        for (int i = 0; i < seqFiles.length; i++)
        {
            Indexer indexer = new Indexer(seqFiles[i], store);
            indexer.setPrimaryKeyName("LOCUS");
            indexer.addSecondaryKey("ACCESSION");

            ChangeTable changeTable = new ChangeTable();

            changeTable.setChanger("LOCUS", new ChangeTable.Changer()
                {
                    public Object change(Object value)
                    {
                        String s = (String) value;
                        int i = s.indexOf(" ");

                        if (i < 0)
                            return s;
                        else
                            return s.substring(0, i);
                    }
                });

            ValueChanger changer = new ValueChanger(indexer, changeTable);
            Parser parser = new Parser();

            while(parser.read(indexer.getReader(),
                              LineSplitParser.GENBANK, changer));
        }

        try
        {
            store.commit();
        }
        catch (CommitFailure ne)
        {
            throw new BioException("Failed to commit new index to file",ne);
        }
    }
}
