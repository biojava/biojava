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

package org.biojava.bio.program.abi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.chromatogram.AbstractChromatogram;
import org.biojava.bio.chromatogram.Chromatogram;
import org.biojava.bio.chromatogram.UnsupportedChromatogramFormatException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.SmallMap;


/**
 * An implementation of {@link org.biojava.bio.chromatogram.Chromatogram} to
 * encapulsulate chromatogram data extracted from the files produced by ABI
 * sequencers, such as the the 377 and the 3700.  The format was described by
 * Clark Tibbetts in his paper "Raw Data File Formats, and the Digital and
 * Analog Raw Data Streams of the ABI PRISM 377 DNA Sequencer."  Available
 * online <kbd><a href="http://www-2.cs.cmu.edu/afs/cs/project/genome/WWW/Papers/clark.html">
 * http://www-2.cs.cmu.edu/afs/cs/project/genome/WWW/Papers/clark.html</a></kbd>
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Richard Holland
 * @see ABIFParser
 */
public class ABIFChromatogram extends AbstractChromatogram implements Serializable {
    public ABIFChromatogram() {
        super();
    }

    /** Create a new ABIF object from a file.
     *  <p>
     *  This method is more efficent than {@link #create(InputStream)}.
     *  </p>
     */
    public static ABIFChromatogram create(File f)
    throws IOException, UnsupportedChromatogramFormatException {
        ABIFChromatogram newOne = new ABIFChromatogram();
        newOne.load(f);
        return newOne;
    }

    /**
     * Create a new ABIF object from a stream of bytes.
     * <p>
     * Due to the non-single-pass design of the ABI format, this method will
     * wrap the InputStream in an {@link org.biojava.utils.io.CachingInputStream}.
     * For this reason, {@link #create(File)} should be preferred.
     * </p>
     * @param in the stream from which to read
     * @return a new ABIFChromatogram object
     * @throws IOException if there is a problem with the underlying stream
     */
    public static ABIFChromatogram create(InputStream in)
    throws IOException, UnsupportedChromatogramFormatException {
        ABIFChromatogram newOne = new ABIFChromatogram();
        newOne.load(in);
        return newOne;
    }

    protected ABIFChromatogram load(File f)
    throws IOException, UnsupportedChromatogramFormatException {
        new Parser(f);
        return this;
    }

    protected ABIFChromatogram load(InputStream in)
    throws IOException, UnsupportedChromatogramFormatException {
        new Parser(in);
        return this;
    }

    protected AbstractChromatogram reverseComplementInstance() {
        return new ABIFChromatogram();
    }

    /**
     * An extension of {@link ABIFParser} that reads the particular fields from
     * the ABIF that contain the chromatogram data and initializes the fields
     * in its enclosing <code>ABIFChromatogram</code> instance.
     */
    protected class Parser extends ABIFParser {
        public Parser(InputStream in)
        throws IOException, UnsupportedChromatogramFormatException {
            super(in);
            parse();
        }

        public Parser(File f)
        throws IOException, UnsupportedChromatogramFormatException {
            super(f);
            parse();
        }

        private final void parse()
        throws IOException, UnsupportedChromatogramFormatException {
            // read filter-wheel-order tag
            char[] fwo_ = new char[4];
            ABIFParser.TaggedDataRecord fwoRec = getDataRecord("FWO_", 1);
            if (fwoRec == null)
                throw new UnsupportedChromatogramFormatException("No FWO_ (1) record in ABIF file, therefore no trace data");
            fwo_[0] = (char) ( (fwoRec.dataRecord >>> 24) & 0xff );
            fwo_[1] = (char) ( (fwoRec.dataRecord >>> 16) & 0xff );
            fwo_[2] = (char) ( (fwoRec.dataRecord >>> 8 ) & 0xff );
            fwo_[3] = (char) ( (fwoRec.dataRecord       ) & 0xff );

            Symbol sym;
            clearTraces();
            for (int i = 0 ; i < 4 ; i++) {
                
                try {
                    sym = ABIFParser.decodeDNAToken(fwo_[i]);
                } catch (IllegalSymbolException ise) {
                    throw new UnsupportedChromatogramFormatException("An unexpected character (" + fwo_[i] +") was found in the FWO_ tag.  Parsing cannot continue.");
                }
                if (!(sym instanceof AtomicSymbol)) {
                    throw new UnsupportedChromatogramFormatException("An unexpected character (" + fwo_[i] +") was found in the FWO_ tag.  Parsing cannot continue.");
                }
                parseTrace((AtomicSymbol) sym, i+9);
            }
            parseBaseCalls();
            getDataAccess().finishedReading();
        }

        private void parseTrace(AtomicSymbol sym, int whichData) throws IOException, UnsupportedChromatogramFormatException {
            TaggedDataRecord dataPtr = getDataRecord("DATA", whichData);
            if (dataPtr.numberOfElements > Integer.MAX_VALUE)
                throw new UnsupportedChromatogramFormatException("Chromatogram has more than " + Integer.MAX_VALUE + " trace samples -- can't handle it");
            int count = (int) dataPtr.numberOfElements;
            int[] trace = new int[count];
            int max = -1;
            setBits(8*dataPtr.elementLength);
            if (dataPtr.elementLength == 2) {
                byte[] shortArray = dataPtr.offsetData;
                int i = 0;
                for (int s = 0; s < shortArray.length; s += 2) {
                    trace[i] =  ((short)((shortArray[s] << 8) | (shortArray[s + 1] & 0xff))) & 0xffff;
                    max = Math.max(trace[i++], max);
                }
            }
            else if (dataPtr.elementLength == 1) {
                byte[] byteArray = dataPtr.offsetData;
                for (int i = 0; i < byteArray.length; i++) {
                    trace[i] = byteArray[i] & 0xff;
                    max = Math.max(trace[i], max);
                }
            }
            else {
                throw new UnsupportedChromatogramFormatException("Only 8- and 16-bit trace samples are supported");
            }
            
            try {
                setTrace(sym, trace, max);
            } catch (IllegalSymbolException ise) {
                throw new BioError("Can't happen", ise);
            }
        }

        private void parseBaseCalls() throws IOException, UnsupportedChromatogramFormatException {
            // do offsets, then call letters
            // offsets are in PLOC1 (we'll use the possibly-edited stream)
            TaggedDataRecord offsetsPtr = getDataRecord("PLOC", 1);
            // call letters are int PBAS1
            TaggedDataRecord basesPtr = getDataRecord("PBAS", 1);
            // these should be equal, but just in case...
            if (offsetsPtr.numberOfElements != basesPtr.numberOfElements)
                throw new BioError("PLOC and PBAS are different lengths.  Can't proceed.");
            if (offsetsPtr.numberOfElements > Integer.MAX_VALUE)
                throw new UnsupportedChromatogramFormatException("Chromatogram has more than " + Integer.MAX_VALUE + " base calls -- can't handle it");
            int count = (int) offsetsPtr.numberOfElements;
            // the list of called bases
            List dna = new ArrayList(count);
            // the list of offsets
            List offsets = new ArrayList(count);
            // start reading offsets, creating SimpleBaseCalls along the way
            if (offsetsPtr.elementLength == 2) {
                byte[] shortArray = offsetsPtr.offsetData;
                IntegerAlphabet integerAlphabet = IntegerAlphabet.getInstance();
                for (int s = 0; s < shortArray.length; s += 2) {
                    offsets.add(integerAlphabet.getSymbol(((short)((shortArray[s] << 8) | (shortArray[s + 1] & 0xff))) & 0xffff));
                }
            }
            else if (offsetsPtr.elementLength == 1) {
                byte[] byteArray = offsetsPtr.offsetData;
                IntegerAlphabet integerAlphabet = IntegerAlphabet.getInstance();
                for (int i = 0 ; i < byteArray.length; i++) {
                    offsets.add(integerAlphabet.getSymbol(byteArray[i] & 0xff));
                }
            }
            else {
                throw new IllegalStateException("Only 8- and 16-bit trace samples are supported");
            }

            // then read the base calls
            try {
                byte[] byteArray = basesPtr.offsetData;
                for (int i = 0; i < byteArray.length; i++) {
                    dna.add(ABIFParser.decodeDNAToken((char) byteArray[i]));
                }
            } catch (IllegalSymbolException ise) {
                throw new BioError("Can't happen", ise);
            }
            // create the base call alignment and set it
            try {
                Map baseCalls = new SmallMap(2);
                baseCalls.put(Chromatogram.DNA, createImmutableSymbolList(DNATools.getDNA(), dna));
                baseCalls.put(Chromatogram.OFFSETS, createImmutableSymbolList(IntegerAlphabet.getInstance(), offsets));
                setBaseCallAlignment(createImmutableAlignment(baseCalls));
            } catch (IllegalAlphabetException iae) {
                throw new BioError("Can't happen", iae);
            } catch (IllegalSymbolException ise) {
                throw new BioError("Can't happen", ise);
            }
        }
    }
}
