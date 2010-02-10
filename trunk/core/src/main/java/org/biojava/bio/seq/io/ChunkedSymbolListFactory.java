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

package org.biojava.bio.seq.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.ChunkedSymbolList;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SimpleSymbolListFactory;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListFactory;
import org.biojava.bio.symbol.SymbolListViews;

/**
 * class that makes ChunkedSymbolLists with the chunks
 * implemented as SymbolLists themselves.
 * <p>
 * The advantage is that those SymbolLists can be packed
 * implementations.
 * <p>
 * You can build a SequenceBuilderFactory to create a packed chunked sequence from
 * an input file without making an intermediate symbol list with:-
 * <pre>
 * public class PackedChunkedListFactory implements SequenceBuilderFactory
 * {
 *   public SequenceBuilder makeSequenceBuilder()
 *   {
 *     return new SequenceBuilderBase() {
 *       private ChunkedSymbolListFactory chunker = new ChunkedSymbolListFactory(new PackedSymbolListFactory(true));
 *
 *       // deal with symbols
 *       public void addSymbols(Alphabet alpha, Symbol[] syms, int pos, int len)
 *         throws IllegalAlphabetException
 *       {
 *         chunker.addSymbols(alpha, syms, pos, len);
 *       }
 *
 *       // make the sequence
 *       public Sequence makeSequence()
 *       {
 *         try {
 *           // make the SymbolList
 *           SymbolList symbols = chunker.makeSymbolList();
 *           seq = new SimpleSequence(symbols, uri, name, annotation);
 *
 *           // call superclass method
 *           return super.makeSequence();
 *         }
 *         catch (IllegalAlphabetException iae) {
 *           throw new BioError("couldn't create symbol list");
 *         }
 *       }
 *     };
 *   }
 * }
 * </pre>
 *
 * <p>
 * Then reading in FASTA files can be done with something like:-
 * <p>
 * <pre>
 * SequenceIterator seqI = new StreamReader(br, new FastaFormat(),
 *     DNATools.getDNA().getTokenization("token"),
 *     new PackedChunkedListFactory() );
 * </pre>
 * <p>
 * Blend to suit taste.
 * <p>
 * Alternatively, you can input Symbols to the factory with addSymbols
 * make the sequence eventually with makeSymbolList.
 * <p>
 * <b>NOTE:</b> An improvement has been introduced where an internal
 * default SymbolList factory is used for small sequences.  This
 * implementation allows for faster SymbolList creation and access
 * for small sequences while allowing a more space-efficient
 * implementation to be selected for large sequences.
 * <p>
 * <b>NOTE:</b> This class is inherantly not threadsafe. You should create one
 * instance for each symbol list you wish to manufacture, and then you should
 * throw that instance away.
 *
 * @author David Huen
 */
public class ChunkedSymbolListFactory
{
    /**
     * operating mode
     */
    private final static int AUTO_SELECT = 1;
    private final static int SUPPLIED_FACTORY = 2;
    private final static int CHUNK_SIZE = 1<<14;

    private final SymbolListFactory userSymListFactory;
    private SymbolListFactory currSymListFactory = null;
    private Alphabet alfa;

    // management variables for chunks
    private Symbol [] headChunk;
    private int headChunkPos = 0;
    private List chunkL = new ArrayList();

    private int symCount = 0;
    private int opMode = SUPPLIED_FACTORY;
    private int threshold = 1<<20;

    // interlocks
    // you can only use symbolAt() or make(), not both.
    private boolean canDoMake = true;

  /**
     * @param symListFactory class which produces the SymbolLists that are used
     *        to store the chunked symbols.
     */
    public ChunkedSymbolListFactory(SymbolListFactory symListFactory)
    {
        opMode = SUPPLIED_FACTORY;
        this.userSymListFactory = symListFactory;
        currSymListFactory = symListFactory;
    }

    /**
     * @param userSymListFactory User-supplied class which produces the SymbolLists
     * that are used to store the chunked symbols (only used when the chunked list
     * to be created is larger than threshold.
     * @param threshold the size of the SymbolList beyond which the userSymListFactory
     * is used.  Below that, the internal default SymbolList factory is used.
     */
    public ChunkedSymbolListFactory(SymbolListFactory userSymListFactory, int threshold)
    {
        opMode = AUTO_SELECT;

        // set threshold to the user-specified value if that's what's requested
        if (threshold > 0) this.threshold = threshold;

        this.userSymListFactory = userSymListFactory;

        // the default factory is the SimpleSymbolList
        currSymListFactory = new SimpleSymbolListFactory();
    }

    /**
     * tool to construct the SymbolList by adding Symbols.
     * Note that this class is not thread-safe.  Also, it
     * can only assemble one SymbolList at a time.  And the
     * composite formed by adding Symbols must not have interstitial
     * missing Symbols.
     */
    public void addSymbols(Alphabet alfa, Symbol[] syms, int pos, int len)
        throws IllegalArgumentException, IllegalAlphabetException
    {
        // lock out make()
        canDoMake = false;

        // check alphabet first
        if (this.alfa == null) {
            this.alfa = alfa;
        }
        else {
            if (this.alfa != alfa) {
                throw new IllegalAlphabetException("Alphabet changed!");
            }
        }

        // increment count
        symCount += len;

        // if count reaches threshold, initiate conversion but do it
        // once only and only if we are in AUTO_SELECT mode
        if ((opMode == AUTO_SELECT)
            && (currSymListFactory != userSymListFactory)
            && (symCount > threshold)) useSuppliedSymListFactory();

        // see if I need to create the first chunk
        if (headChunk == null) {
            headChunk = new Symbol[CHUNK_SIZE];
            headChunkPos = 0;
        }

        // transfer over new symbols
        int ipos = 0;

        while (ipos < len) {
            // is the current chunk full?
            if (headChunkPos == CHUNK_SIZE) {
                // yes.  Verify there are no surprises!
                for (int i= 0; i < CHUNK_SIZE; i++) {
                    if (headChunk[i] == null) {
                        // there's an interstitial null Symbol!
                        throw new IllegalArgumentException("symbols supplied are not tiling contiguously.");
                    }
                }

                // the chunk is alright, lets stash it.
                chunkL.add(currSymListFactory.makeSymbolList(headChunk, CHUNK_SIZE, alfa));
                headChunkPos = 0;
                headChunk = new Symbol[CHUNK_SIZE];
            }

            int read = Math.min(len - ipos, CHUNK_SIZE - headChunkPos);
            System.arraycopy(syms, pos + ipos, headChunk, headChunkPos, read);
            ipos += read;
            headChunkPos += read;
        }
    }

    private void clearState()
    {
        this.canDoMake = true;
        headChunk = null;
        headChunkPos = 0;
        symCount = 0;

        chunkL = new ArrayList();
        alfa = null;

        // if auto-select return to default symbol list factory
        if (opMode == AUTO_SELECT) {
            // the default factory is the SimpleSymbolList
            currSymListFactory = new SimpleSymbolListFactory();
        }
    }

    /**
     * Call this to convert from default SymbolList implementation
     * to user-supplied implementation.
     */
    public void useSuppliedSymListFactory()
    {
        // set the symbolListFactory to the user-supplied one
        currSymListFactory = userSymListFactory;

        // go thru' converting all accumulated chunks to
        // user supplied implementation
        int size = chunkL.size();
        ArrayList temp = new ArrayList(size);

        for (int i = 0; i < size; i++) {
            // get sequence length then transfer contents
            // into a Symbol[] for reencoding.
            SymbolList symList = (SymbolList) chunkL.get(i);
            int symListLen = symList.length();
            Symbol [] symArray;

            if (symList instanceof org.biojava.bio.symbol.SimpleSymbolList) {
                symArray = ((SimpleSymbolList) symList).getSymbolArray();
            }
            else {
                symArray = new Symbol[symListLen];
                for (int j =1; j <= symList.length(); j++) {
                    symArray[j-1] = symList.symbolAt(j);
                }
            }

            // reencode Symbol array: symListLen and CHUNK_SIZE should match except for final chunk.
            // put it into new List
            try {
                temp.add(currSymListFactory.makeSymbolList(symArray, symListLen, alfa));
            }
            catch (IllegalAlphabetException iae) {
                // this should be impossible!!!!!
            }
        }

        // swap over to the new symbol list array
        chunkL = temp;
    }

    /**
     * Converts accumulated Symbols to a SymbolList
     */
    public SymbolList makeSymbolList()
        throws IllegalAlphabetException
    {
        try {
        // we have now acquired all symbols
        // the situation can be;-
        // i) there are no Symbols; then chunkL.size() == 0.
        // ii) one chunk worth, not full: chunkL.size() == 0. headChunkPos != 0;
        // iii) one chunk worth, full: chunkL.size() == 0. headChunkPos == CHUNK_SIZE;
        // iv) multiple chunks, all full: there'll be an unpacked chunk at headChunk;
        // v) multiple chunks, last part filled: chunkL.size() != 0. Unpacked chunk at headChunk;

        // so basically, if chunkL.size() != 0, there will be an unpacked chunk to pack
        // unless headChunkPos == 0;

        // do we have any chunks stashed?
        if (chunkL.size() == 0) {
            // no.  The only chunk if any is at headChunk.
            if (headChunkPos == 0) {
                // no symbols!
                return SymbolListViews.emptyList(alfa);
            }
            else {
                // we do have ONE chunk to deal with.
                // reduce its size as necessary.
                if (headChunkPos < CHUNK_SIZE) {
                    Symbol[] oldChunk = headChunk;
                    headChunk = new Symbol[headChunkPos];
                    System.arraycopy(oldChunk, 0, headChunk, 0, headChunkPos);
                }

                // now return a SubArraySymbolList
                return new SimpleSymbolList(headChunk, headChunkPos, alfa);
            }
        }
        else {
            // we have multiple chunks.

            // do we have an unstashed chunk?
            if (headChunkPos != 0) {
                // yes, let's stash it
                chunkL.add(currSymListFactory.makeSymbolList(headChunk, headChunkPos, alfa));
            }

            // everything we want to put away is now in chunkL
            // create an array to stash the SymbolLists into
            SymbolList [] symListArray = new SymbolList[chunkL.size()];

            // and fill it with the contents of our List
            for (int cnum=0; cnum < chunkL.size(); ++cnum) {
                symListArray[cnum] = (SymbolList) chunkL.get(cnum);
            }

            // now let's return a SymbolList to the user
            int length = (chunkL.size() - 1) * CHUNK_SIZE + headChunkPos;
            return new ChunkedSymbolList(symListArray, CHUNK_SIZE, length, alfa);
        }
        }
        finally {
            clearState();
        }
    }

    /**
     * Method to create a Sequence with a SymbolReader. (does anyone use this???>
     *
     * readme: As of 12/11/03, there are no references to this method in the codebase
     */
    public SymbolList make(SymbolReader sr)
        throws IOException, IllegalSymbolException, IllegalAlphabetException, BioException
    {
        // check interlock
        if(!canDoMake) throw new BioException("you can't use make() and addSymbol() simultaneously.");

        chunkL = new ArrayList();
        headChunk = new Symbol[CHUNK_SIZE];
        headChunkPos = 0;
        alfa = sr.getAlphabet();

        while (sr.hasMoreSymbols()) {
            // is chunk full?
            if (headChunkPos == CHUNK_SIZE) {
                // yes.  Stash away this chunk in packed form.
                chunkL.add(currSymListFactory.makeSymbolList(headChunk, CHUNK_SIZE, sr.getAlphabet()));
                headChunkPos = 0;
            }

            // grab up all available Symbols up to end of chunk.
            int read = sr.readSymbols(headChunk, headChunkPos, CHUNK_SIZE - headChunkPos);
            headChunkPos += read;
        }

        clearState();
        return makeSymbolList();
    }
}
