/*
 * DummyRichSequenceHandler.java
 *
 * Created on March 7, 2006, 3:12 PM
 */

package org.biojavax.bio.seq;

import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.io.ChunkedSymbolListFactory;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.PackedSymbolListFactory;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichLocation.Tools;

/**
 *
 * @author Richard Holland
 * @since 1.5
 */
public class DummyRichSequenceHandler implements RichSequenceHandler {
    /**
     * {@inheritDoc}
     */
    public void edit(RichSequence seq, Edit edit) throws IndexOutOfBoundsException, IllegalAlphabetException, ChangeVetoException {
        seq.getInternalSymbolList().edit(edit);
    }
    
    /**
     * {@inheritDoc}
     */
    public Symbol symbolAt(RichSequence seq, int index) throws IndexOutOfBoundsException {
        if (seq.getCircular()) index = RichLocation.Tools.modulateCircularIndex(index,seq.length());
        return seq.getInternalSymbolList().symbolAt(index);
    }
    
    /**
     * {@inheritDoc}
     */
    public List toList(RichSequence seq) { return seq.getInternalSymbolList().toList();}
    
    /**
     * {@inheritDoc}
     */
    public String subStr(RichSequence seq, int start, int end) throws IndexOutOfBoundsException {
        return seq.getInternalSymbolList().subList(start, end).seqString();
    }
    
    /**
     * {@inheritDoc}
     */
    public SymbolList subList(RichSequence seq, int start, int end) throws IndexOutOfBoundsException {
        if (seq.getCircular()) {
            try {
                int[] modLocation = Tools.modulateCircularLocation(start,end,seq.length());
                int modStart = modLocation[0];
                int modEnd = modLocation[1];
                int modLength = (modEnd - modStart)+1;
                int seqLength = seq.length();
                if (modStart==0) modStart = seqLength;
                if (modEnd==0) modEnd = seqLength;
                // Use the packed symbol factory
                ChunkedSymbolListFactory symsf = new ChunkedSymbolListFactory(new PackedSymbolListFactory());
                if (modEnd>seqLength) {
                    // add it in chunks
                    int remaining = modLength;
                    int chunkSize = (seqLength-modStart)+1;
                    //   add modStart -> seqLength
                    symsf.addSymbols(
                            seq.getAlphabet(),
                            (Symbol[])seq.getInternalSymbolList().subList(modStart,seqLength).toList().toArray(new Symbol[0]),
                            0,
                            chunkSize);
                    remaining -= chunkSize;
                    //   repeat add seqLength
                    while (remaining > seqLength) {
                        chunkSize = seqLength;
                        symsf.addSymbols(
                                seq.getAlphabet(),
                                (Symbol[])seq.getInternalSymbolList().subList(1,seqLength).toList().toArray(new Symbol[0]),
                                0,
                                chunkSize);
                        remaining -= chunkSize;
                    }
                    //   add 0 -> remaining
                    chunkSize = remaining;
                    symsf.addSymbols(
                            seq.getAlphabet(),
                            (Symbol[])seq.getInternalSymbolList().subList(1,chunkSize).toList().toArray(new Symbol[0]),
                            0,
                            chunkSize);
                } else {
                    //   add modStart->modEnd
                    symsf.addSymbols(
                            seq.getAlphabet(),
                            (Symbol[])seq.getInternalSymbolList().subList(modStart,modEnd).toList().toArray(new Symbol[0]),
                            0,
                            modLength);
                }
                return symsf.makeSymbolList();
            } catch (IllegalAlphabetException e) {
                throw new RuntimeException("Don't understand our own alphabet?",e);
            }
        } else return seq.getInternalSymbolList().subList(start, end);
    }
    
    /**
     * {@inheritDoc}
     */
    public String seqString(RichSequence seq) { return seq.getInternalSymbolList().seqString(); }
        
    /**
     * {@inheritDoc}
     */
    public Iterator iterator(RichSequence seq) {return seq.getInternalSymbolList().iterator(); }
}
