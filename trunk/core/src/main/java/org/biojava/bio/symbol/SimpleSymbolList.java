/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public License.  This should
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


package org.biojava.bio.symbol;

import java.io.Serializable;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.io.SeqIOAdapter;
import org.biojava.bio.seq.io.StreamParser;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * Basic implementation of SymbolList.  This
 * is currently backed by a normal Java array.
 *  <p>
 * SimpleSymbolList is now editable. edit() has been implemented
 * in a way that edits are relatively inefficient, but symbolAt() is
 * very efficient.
 * </p>
 * <p>
 * A new constructor SimpleSymbolList(SymbolParser,String) has
 * been added so you can now simply turn a String into a SymbolList.
 * This is mostly to provide a simple way to create a SymbolList for
 * people just trying to get their feet wet. So here is an example.
 * </p>
 * <code>
 * String seqString = "gaattc";
 * FiniteAlphabet dna = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
 * SymbolParser parser = dna.getTokenization("token");
 * SymbolList mySl = new SimpleSymbolList (parser,seqString);
 * System.out.println("Look at my sequence " + mySl.seqString());
 * </code>
 * <p>
 * with the right parser you should be able to make a protein sequence
 * from the String "AspAlaValIleAsp"
 * </p>
 * <p>
 * subList() is implemented such that subLists are views of the original until
 * such time as the underlying SymbolList is edited in a way that would modify
 * the subList, at which point the subList gets its own array of Symbols and
 * does not reflect the edit to the original. When subList() is called on another
 * subList (which is a veiw SimpleSymbolList) the new SimpleSymbolList is a view
 * of the original, not the subList.
 * </p>
 *
 * @author Thomas Down
 * @author David Waring
 * @author David Huen (another constructor)
 * @author George Waldon
 */

public class SimpleSymbolList extends AbstractSymbolList implements ChangeListener, Serializable {
    private static final long serialVersionUID = -9015317520644706924L;
        
    private static int instanceCount;

    private static final int INCREMENT = 100;

    private Alphabet alphabet;
    private Symbol[] symbols;
    private int length;
    private boolean isView;  // this is for subList which returns a view onto a SimpleSymbolList until either is edited
    private int viewOffset;  // offset of the veiw to the original
    private SymbolList referenceSymbolList; // the original SymbolList subLists of views become sublists of original

    private void addListener() {
        incCount();
        alphabet.addChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);
    }

    private static synchronized int incCount() {
        return ++instanceCount;
    }

    protected void finalize() throws Throwable {
        super.finalize();
        // System.err.println("Finalizing a SimpleSymbolList: " + decCount());
        alphabet.removeChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);
        if (isView){
            referenceSymbolList.removeChangeListener(this);
        }
    }

    /**
     * Construct an empty SimpleSymbolList.
     *
     * @param alpha The alphabet of legal symbols in this list.
     */

    public SimpleSymbolList(Alphabet alpha) {
        this.alphabet = alpha;
        this.length = 0;
        this.symbols = new Symbol[INCREMENT];
        this.isView = false;
        this.viewOffset = 0;
        addListener();
    }

    /**
     * Construct a SymbolList containing the symbols in the specified list.
     *
     * @param alpha The alphabet of legal symbols for this list.
     * @param rList A Java List of symbols.
     *
     * @throws IllegalSymbolException if a Symbol is not in the specified alphabet.
     * @throws ClassCastException if rList contains objects which do not implement Symbol.
     */

    public SimpleSymbolList(Alphabet alpha, List rList)
        throws IllegalSymbolException
    {
        this.alphabet = alpha;
        this.length = rList.size();
        symbols = new Symbol[length];
        int pos = 0;
        for (Iterator i = rList.iterator(); i.hasNext(); ) {
            symbols[pos] = (Symbol) i.next();
            alphabet.validate(symbols[pos]);
            pos++;
        }
        this.isView = false;
        this.viewOffset = 0;
        addListener();
    }

    /**
     * Construct a SymbolList from a string.
     *
     * @param parser A SymbolParser for whatever your string is -- e.g. alphabet.getParser("token").
     * @param seqString A Java List of symbols.
     *
     * @throws IllegalSymbolException if a Symbol is not in the specified alphabet.
     */

    public SimpleSymbolList(SymbolTokenization parser, String seqString)
        throws IllegalSymbolException
    {
        if (parser.getTokenType() == SymbolTokenization.CHARACTER) {
            symbols = new Symbol[seqString.length()];
        } else {
            symbols = new Symbol[INCREMENT];
        }
        char[] charArray = new char[1024];
        int segLength = seqString.length();
        StreamParser stParser = parser.parseStream(new SSLIOListener());
        int charCount = 0;
        int chunkLength;
        while (charCount < segLength) {
            chunkLength = Math.min(charArray.length, segLength - charCount);
            seqString.getChars(charCount, charCount + chunkLength, charArray, 0);
            stParser.characters(charArray, 0, chunkLength);
            charCount += chunkLength;
        }
        stParser.close();

        this.alphabet = parser.getAlphabet();
        this.isView = false;
        this.viewOffset = 0;
        addListener();
    }

    /**
     * Construct a copy of an existing SymbolList.
     *
     * @param sl the list to copy.
     */

    public SimpleSymbolList(SymbolList sl) {
        this.alphabet = sl.getAlphabet();
        this.length = sl.length();
        symbols = new Symbol[length];
        for (int i = 0; i < length; ++i) {
            symbols[i] = sl.symbolAt(i + 1);
        }
        this.isView = false;
        this.viewOffset = 0;
        addListener();
    }

    /**
     * Construct a SimpleSymbolList given the Symbol array that backs it.
     * Used primarily with the chunked SymbolList builder but could be used
     * elsewhere too.
     */
    public SimpleSymbolList(Symbol [] symbols, int length, Alphabet alphabet)
    {
        this.symbols = symbols;
        this.length = length;
        this.alphabet = alphabet;
        this.isView = false;
        this.viewOffset = 0;
        addListener();
    }

    /**
     * Construct construct a SimpleSymbolList that is a veiw of the original.
     *    this is used by subList();
     *
     * @param orig -- the original SimpleSymbolList that this is a view of.
     * @param start -- first base in new SymbolList
     * @param end -- last base in new SymbolList
     */

    private SimpleSymbolList(SimpleSymbolList orig, int start, int end) {
        this.alphabet = orig.alphabet;
        this.symbols = orig.symbols;
        this.length = end - start + 1;
        this.isView = true;
        this.viewOffset = start -1;
        this.referenceSymbolList = orig;
        addListener();
    }

    /**
     * Get the alphabet of this SymbolList.
     */

    public Alphabet getAlphabet() {
      return alphabet;
    }

    /**
     * Get the length of this SymbolList.
     */

    public int length() {
      return length;
    }

    /**
     * Find a symbol at a specified offset in the SymbolList.
     *
     * @param pos Position in biological coordinates (1..length)
     */

    public Symbol symbolAt(int pos) {
//        if (pos > length || pos < 1) {
//          throw new IndexOutOfBoundsException(
//            "Can't access " + pos +
//            " as it is not within 1.." + length
//          );
//        }
      // fixme: I have added this check back in a different way as the index
      // system flips from arrays to symbols - we need this detailed debug
      // messaging - anybody want to performance check with/without the try?
      try {
        return symbols[viewOffset + pos - 1];
      } catch (IndexOutOfBoundsException e) {
        throw new IndexOutOfBoundsException(
                "Index must be within [1.." + length() + "] : " + pos);
      }
    }


    /**
    * create a subList of the original, this will be a view until
    * either the original symbolList or the sublist is edited
    */

    public SymbolList subList(int start, int end){
        if (start < 1 || end > length()) {
            throw new IndexOutOfBoundsException(
                      "Sublist index out of bounds " + length() + ":" + start + "," + end
                      );
        }

        if (end < start) {
            throw new IllegalArgumentException(
                "end must not be lower than start: start=" + start + ", end=" + end
                );
        }

        SimpleSymbolList sl = new SimpleSymbolList(this,viewOffset+start,viewOffset+end);
        if (isView){
            referenceSymbolList.addChangeListener(sl);
        }else{
            this.addChangeListener(sl);
        }
        return sl;
    }
    /**
    * Apply and edit to the SymbolList as specified by Edit.
    * <p>
    *   edit() is now supported using the ChangeEvent system. SubLists do NOT reflect edits.
    * </p>
    */

    public synchronized void edit(Edit edit)throws IndexOutOfBoundsException, IllegalAlphabetException,ChangeVetoException {
        ChangeSupport cs;
        ChangeEvent cevt;
        Symbol[] dest;
        int newLength;

        // first make sure that it is in bounds
        if ((edit.pos + edit.length > length +1 ) || (edit.pos <= 0) || edit.length < 0){
            throw new IndexOutOfBoundsException();
        }
        // make sure that the symbolList is of the correct alphabet
        if (( edit.replacement.getAlphabet() != alphabet) &&  (edit.replacement != SymbolList.EMPTY_LIST)){
            throw new IllegalAlphabetException();
        }

        // give the listeners a change to veto this
         // create a new change event ->the EDIT is a static final variable of type ChangeType in SymbolList interface
        cevt = new ChangeEvent(this, SymbolList.EDIT, edit);
        cs = getChangeSupport(SymbolList.EDIT);
        synchronized(cs) {
            // let the listeners know what we want to do
            cs.firePreChangeEvent(cevt);

            // if nobody complained lets continue
            // if we are a view we convert to a real SimpleSymbolList
            if (isView){
                makeReal();
            }
            // now for the edit
            int posRightFragInSourceArray5 = edit.pos + edit.length - 1;
            int rightFragLength = length - posRightFragInSourceArray5;
            int posRightFragInDestArray5 = posRightFragInSourceArray5 + edit.replacement.length() - edit.length;
            int posReplaceFragInDestArray5 = edit.pos - 1;
            int replaceFragLength = edit.replacement.length();

            if ((length + replaceFragLength - edit.length) > symbols.length){
                // extend the array
                dest = new Symbol[(length + replaceFragLength - edit.length + INCREMENT)];

                // copy symbols before the edit no need to do this if we didn't have to build a new array
                System.arraycopy(symbols,0,dest,0,(edit.pos -1));
            }else{
                dest = symbols;  // array copy works when copying from an array to itself
            }

            // copy the symbols after the edit
            if (rightFragLength > 0){
                System.arraycopy(symbols, posRightFragInSourceArray5, dest, posRightFragInDestArray5,rightFragLength);
            }
            // copy the symbols within the edit
            for (int i = 1; i <= replaceFragLength; i++){
                dest[posReplaceFragInDestArray5 + i - 1] = edit.replacement.symbolAt(i);
            }

            // if there was a net deletion we have to get rid of the remaining symbols
            newLength = length + replaceFragLength - edit.length;
            for (int j = newLength; j < length; j++){
                dest[j] = null;
            }
            length = newLength;
            symbols = dest;
            cs.firePostChangeEvent(cevt);
        }
    }

    /**
    *  On preChange() we convert the SymolList to a non-veiw version, giving it its own copy of symbols
    */

    public void preChange(ChangeEvent cev) throws ChangeVetoException{

        // lets not bother making any changes if the edit would not effect us or our children
        Object change = cev.getChange();
        if( (change != null) && (change instanceof Edit) ) {
            Edit e = (Edit)change;
            if (e.pos > (viewOffset + length)){
                return;
            }
            if ((e.pos < viewOffset) && (e.length - e.replacement.length() == 0)){
                return;
            }

        // subLists of views are listeners to the original so we don't have to forward the message
        makeReal();
        }
    }

    // we don't do anything on the postChange we don't want to reflect the changes
    public void postChange(ChangeEvent cev){
    }

    /**
    *  Converts a view symbolList to a real one
    *  that means it gets its own copy of the symbols array
    */
    private void makeReal(){
        if(isView){
            Symbol[] newSymbols = new Symbol[length];
            System.arraycopy (symbols,viewOffset,newSymbols, 0, length);
            this.symbols = newSymbols;
            this.isView = false;
            this.viewOffset = 0;
            referenceSymbolList.removeChangeListener(this);
            referenceSymbolList = null;
        }
    }


    /**
     * Add a new Symbol to the end of this list.
     *
     * @param sym Symbol to add
     * @throws IllegalSymbolException if the Symbol is not in this list's alphabet
     */

      public void addSymbol(Symbol sym)
          throws IllegalSymbolException, ChangeVetoException
      {
          try {
              SymbolList extraSymbol = new SimpleSymbolList(getAlphabet(), Collections.nCopies(1, sym));
              edit(new Edit(length() + 1, 0, extraSymbol));
          } catch (IllegalAlphabetException ex) {
              throw new IllegalSymbolException(ex, sym, "Couldn't add symbol");
          } catch (IndexOutOfBoundsException ex) {
              throw new BioError("Assertion failure: couldn't add symbol at end of list");
          }
      }

    /**
     * Return the Java Symbol[] array that backs this object.
     * primarily used to accelerate reconstruction of symbol lists
     * in the packed chunked symbol list implementation.
     */
    public Symbol [] getSymbolArray()
    {
        return symbols;
    }

    /**
     * Simple inner class for channelling sequence notifications from
     * a StreamParser.
     */

    private class SSLIOListener extends SeqIOAdapter {
        public void addSymbols(Alphabet alpha,Symbol[] syms,int start, int length){
            if(symbols.length < SimpleSymbolList.this.length + length) {
                Symbol[] dest;
                dest = new Symbol [((int) (1.5 * SimpleSymbolList.this.length)) + length];
                System.arraycopy(symbols, 0, dest, 0, SimpleSymbolList.this.length);
                System.arraycopy(syms, start, dest, SimpleSymbolList.this.length, length);
                symbols = dest;
            }else{
                System.arraycopy(syms, start, symbols, SimpleSymbolList.this.length, length);
            }

            SimpleSymbolList.this.length += length;
        }
    }
}
