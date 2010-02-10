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

package org.biojava.bio.program.ssaha;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Symbol;

public interface SequenceStreamer {
    public boolean hasNext();
    public void streamNext(SeqIOListener listener) throws IOException, BioException;
    public void reset() throws BioException;

    public static class SequenceDBStreamer implements SequenceStreamer {
	private SequenceDB seqDB;
	private SequenceIterator si;

	public SequenceDBStreamer(SequenceDB seqDB) {
	    this.seqDB = seqDB;
	    this.si = seqDB.sequenceIterator();
	}

	public boolean hasNext() {
	    return si.hasNext();
	}

	public void reset() {
	    si = seqDB.sequenceIterator();
	}

	public void streamNext(SeqIOListener listener)
	    throws BioException
	{
	    Sequence seq = si.nextSequence();
	    System.err.println("Streaming " + seq.getName());

	    listener.startSequence();
	    listener.setName(seq.getName());
	    listener.setURI(seq.getURN());
	    Symbol[] syms = new Symbol[4096];
	    int pos = 1;
	    int spos = 0;
	    while (pos <= seq.length()) {
		syms[spos++] = seq.symbolAt(pos++);
		if (spos == syms.length || pos > seq.length()) {
		    listener.addSymbols(seq.getAlphabet(), syms, 0, spos);
		    spos = 0;
		}
	    }
	    listener.endSequence();
	}
    }

    public static class FileStreamer implements SequenceStreamer {
	private final List fileList;
	private final SequenceFormat format;
	private final SymbolTokenization toke;
	private Iterator fileIterator;
	private BufferedReader currentStream = null;

	public FileStreamer(SequenceFormat format, SymbolTokenization toke, List files) {
	    this.format = format;
	    this.fileList = files;
	    this.toke = toke;
	    fileIterator = fileList.iterator();
	}

	public FileStreamer(SequenceFormat format, SymbolTokenization toke, File f) {
	    this(format, toke, Collections.singletonList(f));
	}

	public void reset() {
	    currentStream = null;
	    fileIterator = fileList.iterator();
	}

	public boolean hasNext() {
	    return (currentStream != null || fileIterator.hasNext());
	}

	public void streamNext(SeqIOListener listener)
	    throws BioException, IOException
	{
	    if (currentStream == null) {
		currentStream = new BufferedReader(new FileReader((File) fileIterator.next()));
	    }
	    boolean more = format.readSequence(currentStream, toke, listener);
	    if (!more) {
		currentStream.close();
		currentStream = null;
	    }
	}
    }
}
