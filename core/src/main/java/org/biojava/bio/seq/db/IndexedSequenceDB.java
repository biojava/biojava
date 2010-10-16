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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SequenceBuilder;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SequenceBuilderFilter;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ParseErrorListener;
import org.biojava.utils.ParseErrorSource;
import org.biojava.utils.io.CountedBufferedReader;
import org.biojava.utils.io.RandomAccessReader;

/**
 * <p>
 * This class implements SequenceDB on top of a set of sequence files
 * and sequence offsets within these files.
 * </p>
 *
 * <p> This class is primarily responsible for managing the sequence
 * IO, such as calculating the sequence file offsets, and parsing
 * individual sequences based upon file offsets. The actual persistant
 * storage of all this information is delegated to an instance of
 * <code>IndexStore</code>, such as TabIndexStore.
 * </p>
 *
 * <pre>
 * // create a new index store and populate it
 * // this may take some time
 * TabIndexStore indexStore = new TabIndexStore(
 *   storeFile, indexFile, dbName,
 *   format, sbFactory, symbolParser );
 * IndexedSequenceDB seqDB = new IndexedSequenceDB(indexStore);
 *
 * for(int i = 0; i < files; i++) {
 *   seqDB.addFile(files[i]);
 * }
 *
 * // load an existing index store and fetch a sequence
 * // this should be quite quick
 * TabIndexStore indexStore = TabIndexStore.open(storeFile);
 * SequenceDB seqDB = new IndexedSequenceDB(indexStore);
 * Sequence seq = seqDB.getSequence(id);
 * </pre>
 * 
 * <p>
 * Note: We may be able to improve the indexing speed further by
 * discarding all feature creation & annotation requests during index
 * parsing.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Keith James
 * @see org.biojava.bio.seq.db.TabIndexStore
 */

public final class IndexedSequenceDB extends AbstractSequenceDB
    implements SequenceDB, Serializable
{
    private final IDMaker idMaker;
    private final IndexStore indexStore;

    /**
     * Create an IndexedSequenceDB by specifying both the IDMaker and
     * IndexStore used.
     *
     * <p>
     * The IDMaker will be used to calculate the ID for each
     * Sequence. It will delegate the storage and retrieval of the
     * sequence offsets to the IndexStore.
     *
     * @param idMaker  the IDMaker used to calculate Sequence IDs
     * @param indexStore the IndexStore delegate
     */
    public IndexedSequenceDB(IDMaker idMaker, IndexStore indexStore) {
      this.idMaker = idMaker;
      this.indexStore = indexStore;
    }

    /**
     * Create an IndexedSequenceDB by specifying IndexStore used.
     *
     * <p>
     * IDMaker.byName will be used to calculate the ID for each
     * Sequence. It will delegate the storage and retrieval of the
     * sequence offsets to the IndexStore.
     *
     * @param indexStore the IndexStore delegate
     */
    public IndexedSequenceDB(IndexStore indexStore) {
      this(IDMaker.byName, indexStore);
    }

    /**
     * Retrieve the IndexStore.
     *
     * @return the IndexStore delegate
     */
    public IndexStore getIndexStore() {
      return indexStore;
    }

    /**
     * Add sequences from a file to the sequence database. This method
     * works on an "all or nothing" principle. If it can successfully
     * interpret the entire file, all the sequences will be read
     * in. However, if it encounters any problems, it will abandon the
     * whole file; an IOException will be thrown.  Multiple files may
     * be indexed into a single database. A BioException will be
     * thrown if it has problems understanding the sequences.
     *
     * @param seqFile the file containing the sequence or set of sequences
     * @throws BioException if for any reason the sequences can't be read
     *         correctly
     * @throws ChangeVetoException if there is a listener that vetoes adding
     *         the files
     */

    public void addFile(File seqFile)
        throws IllegalIDException, BioException, ChangeVetoException
    {
      boolean completed = false; // initially assume that we will fail
      try {
        seqFile = seqFile.getAbsoluteFile();
        CountedBufferedReader bReader = new CountedBufferedReader(new FileReader(seqFile));
        SequenceFormat format = indexStore.getFormat();
        SymbolTokenization symParser = indexStore.getSymbolParser();
        SequenceBuilderFactory sbFact = indexStore.getSBFactory();
        long pos = bReader.getFilePointer();
        boolean hasNextSequence = true;
        while(hasNextSequence) {
          SequenceBuilder sb = new ElideSymbolsSequenceBuilder(sbFact.makeSequenceBuilder());
          hasNextSequence = format.readSequence(bReader, symParser, sb);
          Sequence seq = sb.makeSequence();
          String id = idMaker.calcID(seq);
          long oldPos = pos;
          pos = bReader.getFilePointer();
          indexStore.store(new SimpleIndex(seqFile, oldPos, (int) (pos - oldPos), id));
        }

        if(!hasListeners()) {
          indexStore.commit();
        } else {
            ChangeEvent ce = new ChangeEvent(
                this,
                SequenceDB.SEQUENCES
            );
            ChangeSupport changeSupport = getChangeSupport(SequenceDB.SEQUENCES);
            synchronized(changeSupport) {
              changeSupport.firePreChangeEvent(ce);
              indexStore.commit();
              changeSupport.firePostChangeEvent(ce);
            }
        }
        completed = true; // we completed succesfuly
      } catch (IOException ioe) {
        throw new BioException("Failed to read sequence file",ioe);
      } finally {
        if(!completed) { // if there was a failure, discard changes
          indexStore.rollback();
        }
      }
    }

    /**
     * SequenceBuilder implementation that explicitly discards the Symbols
     *
     * @author Thomas Down
     * @author Matthew Pocock
     * @author Mark Schreiber
     */
    private static class ElideSymbolsSequenceBuilder
        extends SequenceBuilderFilter
        implements ParseErrorSource {

        private Vector listeners = new Vector();

        public synchronized void removeParseErrorListener(ParseErrorListener p){
          if(listeners.contains(p))
            listeners.remove(p);
        }

        public synchronized void addParseErrorListener(ParseErrorListener p){
          if(! listeners.contains(p)){
            listeners.add(p);
          }
        }

        public ElideSymbolsSequenceBuilder(SequenceBuilder delegate) {
            super(delegate);
        }

        /**
         * Just ignore the symbols
         */
        public void addSymbols(Alphabet alpha, Symbol[] syms, int start, int length) {
        }
    }

    /**
     * Get the name of this sequence database. The name is retrieved
     * from the IndexStore delegate.
     *
     * @return the name of the sequence database, which may be null.
     */
    public String getName() {
        return indexStore.getName();
    }

    public Sequence getSequence(String id)
    throws IllegalIDException, BioException
    {
        try
        {
            Index indx = indexStore.fetch(id);

            RandomAccessReader rar =
                new RandomAccessReader(new RandomAccessFile(indx.getFile(), "r"));

            long toSkip = indx.getStart();
            if (toSkip > rar.length())
                throw new BioException("Reached end of file");
            rar.seek(toSkip);

            SequenceBuilder sb =
                indexStore.getSBFactory().makeSequenceBuilder();

            indexStore.getFormat().readSequence(new BufferedReader(rar),
                                                indexStore.getSymbolParser(),
                                                sb);
            Sequence seq = sb.makeSequence();

            rar.close();
            return seq;
        }
        catch (IOException ioe)
        {
            throw new BioException("Couldn't grab region of file",ioe);
        }
    }

    public SequenceIterator sequenceIterator() {
        return new SequenceIterator() {
            private Iterator idI = indexStore.getIDs().iterator();

            public boolean hasNext() {
                return idI.hasNext();
            }

            public Sequence nextSequence() throws BioException {
                return getSequence((String) idI.next());
            }
        };
    }

    public Set ids() {
        return indexStore.getIDs();
    }


}
