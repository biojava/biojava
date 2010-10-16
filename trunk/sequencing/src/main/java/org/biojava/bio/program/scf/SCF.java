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

package org.biojava.bio.program.scf;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.TreeMap;

import org.biojava.bio.BioError;
import org.biojava.bio.chromatogram.AbstractChromatogram;
import org.biojava.bio.chromatogram.Chromatogram;
import org.biojava.bio.chromatogram.UnsupportedChromatogramFormatException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.utils.SmallMap;



/**
 * A {@link org.biojava.bio.chromatogram.Chromatogram} as loaded from an
 * SCF v2 or v3 file.  Also loads and exposes the SCF format's "private data"
 * and "comments" sections.  The quality values from the SCF are stored as
 * additional sequences on the base call alignment. The labels are the
 * <code>PROB_</code>* constants in this class.
 * The values are {@link org.biojava.bio.symbol.IntegerAlphabet.IntegerSymbol}
 * objects in the range 0 to 255.
 *
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 */
public class SCF extends AbstractChromatogram {
    private byte[] privateData;
    private Properties comments;
    
    private static final IntegerAlphabet.SubIntegerAlphabet
            PROBABILITY_ALPHABET =
            IntegerAlphabet.getSubAlphabet(0, 255);
    
    /** Represents the maximum unsigned value
     * of a byte for wrapping purposes */
    public static final int BYTE_MAX_VALUE =
            256;
    
    /** Represents the maximum unsigned value
     * of a short for wrapping purposes */
    public static final int SHORT_MAX_VALUE =
            65536;
    
    /** Base call alignment sequence label for the probability that call
     * should be A. */
    public static final String PROB_NUC_A = "quality-a";
    /** Base call alignment sequence label for the probability that call
     * should be C. */
    public static final String PROB_NUC_C = "quality-c";
    /** Base call alignment sequence label for the probability that call
     * should be G. */
    public static final String PROB_NUC_G = "quality-g";
    /** Base call alignment sequence label for the probability that call
     * should be T. */
    public static final String PROB_NUC_T = "quality-t";
    /**
     * Base call alignment sequence label for the substitution
     * probability. In versions of the SCF spec before 3.10, this is called
     * spareQual[0].
     */
    public static final Object PROB_SUBSTITUTION =
            "substitution-probability";
    /**
     * Base call alignment sequence label for the overcall probability.
     * In versions of the SCF spec before 3.10, this is called
     * spareQual[1].
     */
    public static final Object PROB_OVERCALL     = "overcall-probability";
    /**
     * Base call alignment sequence label for the undercall probability.
     * In versions of the SCF spec before 3.10, this is called
     * spareQual[2].
     */
    public static final Object PROB_UNDERCALL    =
            "undercall-probability";
    
    /** Creates a new, completely empty SCF. */
    protected SCF() {
        super();
        comments = new Properties();
    }
    
    public static SCF create(File f)
    throws IOException, UnsupportedChromatogramFormatException {
        SCF out = new SCF();
        out.load(f);
        return out;
    }
    
    public static SCF create(InputStream in, long alreadyRead)
    throws IOException, UnsupportedChromatogramFormatException {
        SCF out = new SCF();
        out.load(in, alreadyRead);
        return out;
    }
    
    protected void load(File f) throws IOException,
            UnsupportedChromatogramFormatException {
        FileInputStream fin = new FileInputStream(f);
        try {
        	load(fin, 0);
        } finally {
        	fin.close();
        }
    }
    
    protected void load(InputStream in, long initOffset)
    throws IOException, UnsupportedChromatogramFormatException {
        SCF.ParserFactory.parse(in, this, initOffset);
    }
    
    /**
     * Returns the comments fields as a {@link Properties} mapping.
     */
    public Properties getComments() { return comments; }
    
    protected AbstractChromatogram reverseComplementInstance() { return
            new SCF(); }
    
    public static IntegerAlphabet.SubIntegerAlphabet
            getProbabilityAlphabet() { return PROBABILITY_ALPHABET; }
    
    /**
     * Overrides {@link
     * AbstractChromatogram#reverseComplementBaseCallList} to
     * support the 7 quality values from the SCF.  These are handled thus:
     * <ul>
     *   <li><code>PROB_SUBSTITUTION</code>, <code>PROB_OVERCALL</code>, and
     *       <code>PROB_UNDERCALL</code> are just reversed &returned.</li>
     *   <li><code>PROB_NUC_</code>* returns the reverse of the quality
     *       sequence for the complement base.</li>
     * </ul>
     */
    protected SymbolList reverseComplementBaseCallList(String label) {
        if (label == PROB_SUBSTITUTION ||
                label == PROB_OVERCALL ||
                label == PROB_UNDERCALL) {
            return
                    SymbolListViews.reverse(this.getBaseCalls().symbolListForLabel(label));
        } else if (label.equals( PROB_NUC_A)) {
            return
                    SymbolListViews.reverse(this.getBaseCalls().symbolListForLabel(PROB_NUC_T));
        } else if (label.equals(PROB_NUC_C)) {
            return
                    SymbolListViews.reverse(this.getBaseCalls().symbolListForLabel(PROB_NUC_G));
        } else if (label.equals(PROB_NUC_G)) {
            return
                    SymbolListViews.reverse(this.getBaseCalls().symbolListForLabel(PROB_NUC_C));
        } else if (label.equals( PROB_NUC_T)) {
            return
                    SymbolListViews.reverse(this.getBaseCalls().symbolListForLabel(PROB_NUC_A));
        } else {
            return super.reverseComplementBaseCallList(label);
        }
    }
    
    /*** INNER CLASSES FOR PARSER ***/
    
    /** Factory class to create the appropriate parser for the given
     * stream.  This decision is based on the version field in the file's header. */
    private static class ParserFactory {
        public static void parse(InputStream in, SCF out, long initOffset)
        throws IOException, UnsupportedChromatogramFormatException {
            DataInputStream din = new DataInputStream(in);
            SCF.Parser parser = createParser(din, out, initOffset);
            parser.parse();
        }
        
        public static SCF.Parser createParser(DataInputStream din, SCF
                out, long initOffset)
                throws UnsupportedChromatogramFormatException, IOException {
            // read the header to find out the version
            long offset = initOffset;
            SCF.Parser.HeaderStruct header =
                    SCF.Parser.HeaderStruct.create(din, initOffset);
            offset = SCF.Parser.HeaderStruct.HEADER_LENGTH;
            out.setBits((int)header.sample_size * 8);
            SCF.Parser parser;
            float version;
            try {
                version = Float.parseFloat(new String(header.version));
            } catch (NumberFormatException e) {
                throw new UnsupportedChromatogramFormatException(
                        "The SCF's version (" + new String(header.version) +
                        ") is not a number");
            }
            if (version < 3.0f && version >= 2.0f) {
                parser = new SCF.V2Parser(din, out, header, offset);
            } else if (version >= 3.0f) {
                parser = new SCF.V3Parser(din, out, header, offset);
            } else {
                throw new UnsupportedChromatogramFormatException(
                        "Only version 2 and version 3 SCFs are supported (not "
                        + new String(header.version));
            }
            return parser;
        }
    }
    
    static interface BaseCallUncertaintyDecoder {
        /**
         * Returns an appropriate Symbol from the DNA alphabet for
         * an encoded byte.
         */
        public Symbol decode(byte call) throws IllegalSymbolException;
    }
    
    /**
     * A BaseCallUncertaintyDecoder that works for type 0 (default) and
     * type 4 (ABI)
     * code sets.
     */
    static class DefaultUncertaintyDecoder implements
            BaseCallUncertaintyDecoder {
        public DefaultUncertaintyDecoder() { }
        public Symbol decode(byte call) throws IllegalSymbolException {
            char c = (char) call;
            switch (c) {
                case 'a': case 'A':
                    return DNATools.a();
                case 'c': case 'C':
                    return DNATools.c();
                case 'g': case 'G':
                    return DNATools.g();
                case 't': case 'T':
                    return DNATools.t();
                case 'n': case 'N':
                    return DNATools.n();
                case 'm': case 'M':
                    return DNATools.m();
                case 'r': case 'R':
                    return DNATools.r();
                case 'w': case 'W':
                    return DNATools.w();
                case 's': case 'S':
                    return DNATools.s();
                case 'y': case 'Y':
                    return DNATools.y();
                case 'k': case 'K':
                    return DNATools.k();
                case 'v': case 'V':
                    return DNATools.v();
                case 'h': case 'H':
                    return DNATools.h();
                case 'd': case 'D':
                    return DNATools.d();
                case 'b': case 'B':
                    return DNATools.b();
                case '-':
                    return DNATools.getDNA().getGapSymbol();
                default:
                    throw new IllegalSymbolException("No Symbol for " +
                            c);
            }
        }
    }
    
    abstract static class Parser {
        protected long offset = 0;
        protected DataInputStream din = null;
        protected HeaderStruct header;
        protected BaseCallUncertaintyDecoder decoder;
        protected SCF out = null;
        protected boolean parsed = false;
        
        Parser(DataInputStream din, SCF out,
                SCF.Parser.HeaderStruct header, long initOffset)
                throws UnsupportedChromatogramFormatException {
            if (din == null)
                throw new IllegalArgumentException(
                        "Can't parse a null inputstream");
            this.din = din;
            
            if (out == null)
                this.out = new SCF();
            else
                this.out = out;
            
            if (header.samples > Integer.MAX_VALUE)
                throw new UnsupportedChromatogramFormatException(
                        "Can't parse an SCF with more than " + Integer.MAX_VALUE + " trace samples");
            
            if (header.bases > Integer.MAX_VALUE)
                throw new UnsupportedChromatogramFormatException(
                        "Can't parse an SCF with more than " + Integer.MAX_VALUE + " called bases");
            
            this.header = header;
            this.decoder = createDecoder(header.code_set);
            this.offset = initOffset;
        }
        
        /**
         * Factory method to create an approriate decoder for the code
         * set. Currently, only a direct interpretation of the encoded byte
         * as an ASCII character is supported (see DefaultUncertaintyDecoder).
         * This interpretation will be used for all code sets, but a
         * warning will be printed if the code set is not known to work with this
         * interpretation.
         */
        private static BaseCallUncertaintyDecoder createDecoder(long
                codeSet) {
            if (codeSet != 0 && codeSet != 4)
                System.err.println("Warning: the code set (" + codeSet +
                        ") is not specifically supported.  (It may still work, though.)");
            return new DefaultUncertaintyDecoder();
        }
        
        public SCF getParsed() {
            if (parsed) return out;
            else        return null;
        }
        
        public void parse() throws IOException,
                UnsupportedChromatogramFormatException {
            parsed = false;
            // sort the sections of the file by ascending offset
            Integer SAMPLES  = new Integer(0),
                    BASES    = new Integer(1),
                    COMMENTS = new Integer(2),
                    PRIVATE  = new Integer(3);
            TreeMap sectionOrder = new TreeMap();
            sectionOrder.put(new Long(header.samples_offset),  SAMPLES);
            sectionOrder.put(new Long(header.bases_offset),    BASES);
            sectionOrder.put(new Long(header.comments_offset), COMMENTS);
            sectionOrder.put(new Long(header.private_offset),  PRIVATE);
            
            for (Iterator it = sectionOrder.keySet().iterator() ;
            it.hasNext() ;) {
                Integer sect = (Integer) sectionOrder.get(it.next());
                if      (sect == SAMPLES)  parseSamples();
                else if (sect == BASES)    parseBases();
                else if (sect == COMMENTS) parseComments();
                else if (sect == PRIVATE)  parsePrivate();
            }
            parsed = true;
        }
        
        protected abstract void parseSamples()  throws IOException,
                UnsupportedChromatogramFormatException;
        protected abstract void parseBases()    throws IOException,
                UnsupportedChromatogramFormatException;
        
        protected void parseComments() throws IOException {
            skipTo(header.comments_offset);
            byte[] raw = new byte[(int)header.comments_size - 1];
            din.read(raw, 0, raw.length);
            BufferedReader r = new BufferedReader(
                    new InputStreamReader(
                    new ByteArrayInputStream(raw),
                    "ISO-8859-1"
                    )
                    );
            String line, key, value;
            int eqIdx = -1;
            while ((line = r.readLine()) != null) {
                eqIdx = line.indexOf('=');
                //added line below to skip any truncated comment fields
                if( eqIdx == -1 ) {
                    continue;
                    
                }
                key = line.substring(0, eqIdx);
                value = line.substring(eqIdx+1);
                out.comments.setProperty(key, value);
            }
        }
        
        protected void parsePrivate() throws IOException {
            if (header.private_size == 0) return;
            skipTo(header.private_offset);
            out.privateData = new byte[(int)header.private_size];
            int privRead = 0;
            int thisRead = 0;
            while (privRead < out.privateData.length && thisRead >= 0) {
                thisRead = din.read(out.privateData, privRead,
                        out.privateData.length - privRead);
               	offset += thisRead;
               	privRead += thisRead;
            }
        }
        
        protected final void skipTo(long newOffset) throws IOException {
            if (newOffset < offset)
                throw new IllegalArgumentException("Can't skip backwards: (newOffset==" + newOffset + ") < (offset==" + offset + ")");
            long skip = newOffset - offset;
            while (skip > 0) {
                int actualSkip = din.skipBytes((int)skip);
                offset += actualSkip;
                skip   -= actualSkip;
            }
        }
        
        /**
         * Does the grunt work of creating the base call alignment from
         * the given lists of bases, offsets, and probabilities.
         * Catches and "Can't happens" all exceptions.
         */
        protected final void createAndSetBaseCallAlignment(List dna, List
                offsets, List[] probs) {
            try {
                Map baseCalls = new SmallMap(9);
                baseCalls.put(Chromatogram.DNA,
                        out.createImmutableSymbolList(DNATools.getDNA(), dna));
                baseCalls.put(Chromatogram.OFFSETS,
                        out.createImmutableSymbolList(IntegerAlphabet.getInstance(), offsets));
                baseCalls.put(PROB_NUC_A,
                        out.createImmutableSymbolList(getProbabilityAlphabet(), probs[0]));
                baseCalls.put(PROB_NUC_C,
                        out.createImmutableSymbolList(getProbabilityAlphabet(), probs[1]));
                baseCalls.put(PROB_NUC_G,
                        out.createImmutableSymbolList(getProbabilityAlphabet(), probs[2]));
                baseCalls.put(PROB_NUC_T,
                        out.createImmutableSymbolList(getProbabilityAlphabet(), probs[3]));
                baseCalls.put(PROB_SUBSTITUTION,
                        out.createImmutableSymbolList(getProbabilityAlphabet(), probs[4]));
                baseCalls.put(PROB_OVERCALL,
                        out.createImmutableSymbolList(getProbabilityAlphabet(), probs[5]));
                baseCalls.put(PROB_UNDERCALL,
                        out.createImmutableSymbolList(getProbabilityAlphabet(), probs[6]));
                out.setBaseCallAlignment(out.createImmutableAlignment(baseCalls));
            } catch (IllegalSymbolException ise) {
                throw new BioError(ise,"Can't happen unless the decoder is returning non-DNA symbols");
            } catch (IllegalAlphabetException iae) {
                throw new BioError(iae,"Can't happen");
            }
        }
        
        private static class HeaderStruct {
            public static final int HEADER_LENGTH = 128;
            
            // SCF spec uses unsigned 32-bit ints, so we'll use longs for simplicity
            public long magic_number;
            public long samples;
            public long samples_offset;
            public long bases;
            public long bases_left_clip;
            public long bases_right_clip;
            public long bases_offset;
            public long comments_size;
            public long comments_offset;
            public char[] version;
            public long sample_size;
            public long code_set;
            public long private_size;
            public long private_offset;
            public long[] spare;
            
            private HeaderStruct() {
                version = new char[4];
                spare = new long[18];
            }
            
            public static HeaderStruct create(DataInputStream din, long
                    initOffset) throws IOException {
                HeaderStruct hs = new HeaderStruct();
                if (initOffset > 4) {
                    throw new IllegalStateException("Can't skip more than four bytes and still have enough info to read header");
                } else if (initOffset == 0) {
                    hs.magic_number = 0xFFFFFFFF & din.readInt();
                } else {
                    hs.magic_number = 0;
                    // skip to the four-byte boundary
                    for (int i = 0 ; i < 4 - initOffset ; i++)
                        din.read();
                }
                hs.samples          = 0xFFFFFFFF & din.readInt();
                hs.samples_offset   = 0xFFFFFFFF & din.readInt();
                hs.bases            = 0xFFFFFFFF & din.readInt();
                hs.bases_left_clip  = 0xFFFFFFFF & din.readInt();
                hs.bases_right_clip = 0xFFFFFFFF & din.readInt();
                hs.bases_offset     = 0xFFFFFFFF & din.readInt();
                hs.comments_size    = 0xFFFFFFFF & din.readInt();
                hs.comments_offset  = 0xFFFFFFFF & din.readInt();
                hs.version[0]       = (char) din.readByte();
                hs.version[1]       = (char) din.readByte();
                hs.version[2]       = (char) din.readByte();
                hs.version[3]       = (char) din.readByte();
                hs.sample_size      = 0xFFFFFFFF & din.readInt();
                hs.code_set         = 0xFFFFFFFF & din.readInt();
                hs.private_size     = 0xFFFFFFFF & din.readInt();
                hs.private_offset   = 0xFFFFFFFF & din.readInt();
                for (int i = 0 ; i < hs.spare.length ; i++) {
                    hs.spare[i] = 0xFFFFFFFF & din.readInt();
                }
                return hs;
            }
        }
    }
    
    private static class V3Parser extends Parser {
        V3Parser(DataInputStream din, SCF out,
                SCF.Parser.HeaderStruct header, long initOffset)
                throws IOException, UnsupportedChromatogramFormatException {
            super(din, out, header, initOffset);
        }
        
        protected void parseSamples() throws IOException,
                UnsupportedChromatogramFormatException {
            skipTo(header.samples_offset);
            
            // load values from file
            int count = (int)header.samples;
            int maxValAllowed = 0;
            if(header.sample_size == 1) {
                maxValAllowed = BYTE_MAX_VALUE;
            } else if(header.sample_size == 2) {
                maxValAllowed = SHORT_MAX_VALUE;
            }
            
            int[][] trace = new int[4][count];
            int[] maxVal = new int[] { Integer.MIN_VALUE,
                    Integer.MIN_VALUE,
                    Integer.MIN_VALUE,
                    Integer.MIN_VALUE };
                    for (int n = 0 ; n < 4 ; n++)
                        readSamplesInto(trace[n]);
                    
                    // stored values are delta-delta values; reprocess into actual values
                    // algorithm cribbed from io_lib's delta_samples function in misc_scf.c
                    int[] p_sample1 = new int[] { 0, 0, 0, 0 };
                    int[] p_sample2 = new int[] { 0, 0, 0, 0 };
                    for (int n = 0 ; n < 4 ; n++) {
                        for (int i = 0 ; i < trace[n].length ; i++) {
                            //New version of the code which takes into account what the
                            //underlying value i.e. is it a byte or is it a short
                            //and slightly rejigged for speed and sanity sakes
                            p_sample1[n] += trace[n][i];
                            if(p_sample1[n] >= maxValAllowed) p_sample1[n] = p_sample1[n] -
                                    maxValAllowed;
                            p_sample2[n] = p_sample1[n] + p_sample2[n];
                            if(p_sample2[n] >= maxValAllowed) p_sample2[n] = p_sample2[n] -
                                    maxValAllowed;
                            trace[n][i] = p_sample2[n];
                            
                            maxVal[n] = Math.max(maxVal[n], trace[n][i]);
                        }
                    }
                    // set into output SCF chromat object
                    try {
                        out.setTrace(DNATools.a(), trace[0], maxVal[0]);
                        out.setTrace(DNATools.c(), trace[1], maxVal[1]);
                        out.setTrace(DNATools.g(), trace[2], maxVal[2]);
                        out.setTrace(DNATools.t(), trace[3], maxVal[3]);
                    } catch (IllegalSymbolException ise) {
                        throw new BioError(ise,"Can't happen");
                    }
        }
        
        protected void readSamplesInto(int[] samps) throws IOException {
            if (header.sample_size == 1) {
                for (int i = 0 ; i < samps.length ; i++) {
                    samps[i] = din.readUnsignedByte();
                    offset += 1;
                }
            } else if (header.sample_size == 2) {
                for (int i = 0 ; i < samps.length ; i++) {
                    samps[i] = din.readUnsignedShort();
                    offset += 2;
                }
            }
        }
        
        protected void parseBases() throws IOException,
                UnsupportedChromatogramFormatException {
            skipTo(header.bases_offset);
            int count = (int) header.bases;
            
            List[] probs = new ArrayList[7];
            for (int i = 0 ; i < 7 ; i++) probs[i] = new ArrayList(count);
            List offsets = new ArrayList(count);
            List dna     = new ArrayList(count);
            
            long tmp;
            try {
                for (int i = 0 ; i < count ; i++) {
                    // read the first set of sequential data (the trace peak offsets)
                    tmp = 0xFFFFFFFF & din.readInt();
                    offset += 4;
                    if (tmp > Integer.MAX_VALUE)
                        throw new
                                UnsupportedChromatogramFormatException(
                                "SCF contains a base with peak offset > " + Integer.MAX_VALUE);
                    offsets.add(IntegerAlphabet.getInstance().getSymbol((int) tmp));
                }
                // read sets of probs for A, C, G, T
                for (int j = 0 ; j < 4 ; j++) {
                    for (int i = 0 ; i < count ; i++) {
                        probs[j].add(getProbabilityAlphabet().getSymbol(din.readByte() & 0xFF));
                        offset++;
                    }
                }
                // read bases
                try {
                    for (int i = 0 ; i < count ; i++) {
                        dna.add(decoder.decode(din.readByte()));
                        offset++;
                    }
                } catch (IllegalSymbolException ise) {
                	UnsupportedChromatogramFormatException ue = new UnsupportedChromatogramFormatException("Base call decoding failure");
                	ue.initCause(ise);
                	throw ue;
                }
                // read 'spare' probs
                for (int j = 4 ; j < 7 ; j++) {
                    for (int i = 0 ; i < count ; i++) {
                        probs[j].add(getProbabilityAlphabet().getSymbol(din.readByte() & 0xFF));
                        offset++;
                    }
                }
            } catch (IllegalSymbolException ise) {
                throw new BioError(ise,"Can't happen unless there's a misdefinition of getProbabilityAlphabet() or IntegerAlphabet");
            }
            // create/set base call list
            createAndSetBaseCallAlignment(dna, offsets, probs);
        }
    } // end SCFv3Parser
    
    private static class V2Parser extends Parser {
        V2Parser(DataInputStream din, SCF out,
                SCF.Parser.HeaderStruct header, long initOffset)
                throws IOException, UnsupportedChromatogramFormatException {
            super(din, out, header, initOffset);
        }
        
        protected void parseSamples() throws IOException {
            int count = (int) header.samples;
            int[][] trace = new int[4][count];
            int[] maxVal = new int[] { Integer.MIN_VALUE,
                    Integer.MIN_VALUE,
                    Integer.MIN_VALUE,
                    Integer.MIN_VALUE };
                    
                    if (header.sample_size == 1) {
                        for (int i = 0 ; i < count ; i++) {
                            for (int n = 0 ; n < 4 ; n++) {
                                trace[n][i] = din.readUnsignedByte();
                                maxVal[n] = Math.max(trace[n][i], maxVal[n]);
                                offset++;
                            }
                        }
                    } else if (header.sample_size == 2) {
                        for (int i = 0 ; i < count ; i++) {
                            for (int n = 0 ; n < 4 ; n++) {
                                trace[n][i] = din.readUnsignedShort();
                                maxVal[n] = Math.max(trace[n][i], maxVal[n]);
                                offset += 2;
                            }
                        }
                    }
                    // set into output SCF chromat object
                    try {
                        out.setTrace(DNATools.a(), trace[0], maxVal[0]);
                        out.setTrace(DNATools.c(), trace[1], maxVal[1]);
                        out.setTrace(DNATools.g(), trace[2], maxVal[2]);
                        out.setTrace(DNATools.t(), trace[3], maxVal[3]);
                    } catch (IllegalSymbolException ise) {
                        throw new BioError(ise,"Can't happen");
                    }
        }
        
        protected void parseBases() throws IOException,
                UnsupportedChromatogramFormatException {
            skipTo(header.bases_offset);
            
            int count = (int) header.bases;
            List[] probs = new ArrayList[7];
            for (int i = 0 ; i < probs.length ; i++) probs[i] = new
                    ArrayList(count);
            List dna = new ArrayList(count);
            List offsets = new ArrayList(count);
            long tmp;
            byte[] probTmp = new byte[7];
            for (int i = 0 ; i < count ; i++) {
                // read the peak index
                tmp = 0xFFFFFFFF & din.readInt();
                offset += 4;
                if (tmp > Integer.MAX_VALUE)
                    throw new UnsupportedChromatogramFormatException(
                            "SCF contains a base with peak offset > " + Integer.MAX_VALUE);
                offsets.add(IntegerAlphabet.getInstance().getSymbol((int)
                tmp));
                
                // read the per-base probabilities
                din.read(probTmp, 0, 4);
                offset += 4;
                // read the actual base called
                try {
                    dna.add(decoder.decode(din.readByte()));
                    offset += 1;
                } catch (IllegalSymbolException ise) {
                	UnsupportedChromatogramFormatException ue = new UnsupportedChromatogramFormatException(
                            "Base call decoding failure");
                	ue.initCause(ise);
                	throw ue;
                }
                // read the spare probability fields
                din.read(probTmp, 4, 3);
                offset += 3;
                try {
                    for (int p = 0 ; p < 7 ; p++)
                        probs[p].add(getProbabilityAlphabet().getSymbol(0xFF & probTmp[p]));
                } catch (IllegalSymbolException ise) {
                    throw new BioError(ise,"Can't happen unless getProbabilityAlphabet() has been misdefined.");
                }
            }
            createAndSetBaseCallAlignment(dna, offsets, probs);
        }
    }
}


