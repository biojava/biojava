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

package org.biojava.bio.seq;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.PairDistribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.seq.impl.SimpleGappedSequence;
import org.biojava.bio.seq.impl.SimpleSequenceFactory;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.ReversibleTranslationTable;
import org.biojava.bio.symbol.SimpleReversibleTranslationTable;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.utils.ChangeVetoException;

/**
 * Useful functionality for processing DNA sequences.
 *
 * @author Matthew Pocock
 * @author Keith James (docs)
 * @author Mark Schreiber
 * @author David Huen
 * @author Richard Holland
 */
public final class DNATools {
  private static final ReversibleTranslationTable complementTable;
  static private final FiniteAlphabet dna;
    private static final SymbolTokenization dnaTokens;

  static private final AtomicSymbol a;
  static private final AtomicSymbol g;
  static private final AtomicSymbol c;
  static private final AtomicSymbol t;
  static private final Symbol n;
  static private final Symbol m;
  static private final Symbol r;
  static private final Symbol w;
  static private final Symbol s;
  static private final Symbol y;
  static private final Symbol k;
  static private final Symbol v;
  static private final Symbol h;
  static private final Symbol d;
  static private final Symbol b;
  static private final SimpleReversibleTranslationTable transcriptionTable;


  static private Map<Symbol, Symbol> symbolToComplement;

  static {
    try {
      dna = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
      dnaTokens = dna.getTokenization("token");
      SymbolList syms = new SimpleSymbolList(dnaTokens, "agctnmrwsykvhdb");
      a = (AtomicSymbol) syms.symbolAt(1);
      g = (AtomicSymbol) syms.symbolAt(2);
      c = (AtomicSymbol) syms.symbolAt(3);
      t = (AtomicSymbol) syms.symbolAt(4);
      n = syms.symbolAt(5);
      m = syms.symbolAt(6);
      r = syms.symbolAt(7);
      w = syms.symbolAt(8);
      s = syms.symbolAt(9);
      y = syms.symbolAt(10);
      k = syms.symbolAt(11);
      v = syms.symbolAt(12);
      h = syms.symbolAt(13);
      d = syms.symbolAt(14);
      b = syms.symbolAt(15);

      symbolToComplement = new HashMap<Symbol, Symbol>();

      // add the gap symbol
      Symbol gap = dna.getGapSymbol();
      symbolToComplement.put(gap, gap);

      // add all other ambiguity symbols
      for(Iterator i = AlphabetManager.getAllSymbols(dna).iterator(); i.hasNext();) {
          Symbol as = (Symbol) i.next();
          FiniteAlphabet matches = (FiniteAlphabet) as.getMatches();
          if (matches.size() > 1) {   // We've hit an ambiguous symbol.
              Set<Symbol> l = new HashSet<Symbol>();
              for(Iterator j = matches.iterator(); j.hasNext(); ) {
                  l.add(complement((Symbol) j.next()));
              }
              symbolToComplement.put(as, dna.getAmbiguity(l));
          }
      }


      complementTable = new DNAComplementTranslationTable();
      
      transcriptionTable = new SimpleReversibleTranslationTable(getDNA(), RNATools.getRNA());
      transcriptionTable.setTranslation(a, RNATools.a());
      transcriptionTable.setTranslation(c, RNATools.c());
      transcriptionTable.setTranslation(g, RNATools.g());
      transcriptionTable.setTranslation(t, RNATools.u());
      
    } catch (Throwable th) {
      throw new BioError("Unable to initialize DNATools", th);
    }
  }

  public static AtomicSymbol a() { return a; }
  public static AtomicSymbol g() { return g; }
  public static AtomicSymbol c() { return c; }
  public static AtomicSymbol t() { return t; }
  public static Symbol n() { return n; }
  public static Symbol m() { return m; }
  public static Symbol r() { return r; }
  public static Symbol w() { return w; }
  public static Symbol s() { return s; }
  public static Symbol y() { return y; }
  public static Symbol k() { return k; }
  public static Symbol v() { return v; }
  public static Symbol h() { return h; }
  public static Symbol d() { return d; }
  public static Symbol b() { return b; }

  
  private DNATools() {
  }
  
  /**
   * Return the DNA alphabet.
   *
   * @return a flyweight version of the DNA alphabet
   */
  public static FiniteAlphabet getDNA() {
    return dna;
  }

  /**
   * Gets the (DNA x DNA) Alphabet
   * @return a flyweight version of the (DNA x DNA) alphabet
   */
  public static FiniteAlphabet getDNAxDNA(){
    return (FiniteAlphabet)AlphabetManager.generateCrossProductAlphaFromName("(DNA x DNA)");
  }

  /**
   * Gets the (DNA x DNA x DNA) Alphabet
   * @return a flyweight version of the (DNA x DNA x DNA) alphabet
   */
  public static FiniteAlphabet getCodonAlphabet(){
    return (FiniteAlphabet)AlphabetManager.generateCrossProductAlphaFromName("(DNA x DNA x DNA)");
  }

  /**
   * Return a new DNA <span class="type">SymbolList</span> for
   * <span class="arg">dna</span>.
   *
   * @param dna a <span class="type">String</span> to parse into DNA
   * @return a <span class="type">SymbolList</span> created form
   *         <span class="arg">dna</span>
   * @throws IllegalSymbolException if <span class="arg">dna</span> contains
   *         any non-DNA characters
   */
  public static SymbolList createDNA(String dna)
  throws IllegalSymbolException {
    SymbolTokenization p = null;
    try {
      p = getDNA().getTokenization("token");
    } catch (BioException e) {
      throw new BioError("Something has gone badly wrong with DNA", e);
    }
    return new SimpleSymbolList(p, dna);
  }

  /**
   * Return a new DNA <span class="type">Sequence</span> for
   * <span class="arg">dna</span>.
   *
   * @param dna a <span class="type">String</span> to parse into DNA
   * @param name a <span class="type">String</span> to use as the name
   * @return a <span class="type">Sequence</span> created form
   *         <span class="arg">dna</span>
   * @throws IllegalSymbolException if <span class="arg">dna</span> contains
   *         any non-DNA characters
   */
  public static Sequence createDNASequence(String dna, String name)
  throws IllegalSymbolException {
    //should I be calling createGappedDNASequence?
    if(dna.indexOf('-') != -1 || dna.indexOf('~') != -1){//there is a gap
        return createGappedDNASequence(dna, name);
    }
    // No need to wrap in try-catch as throws IllegalSymbolException already.
    //try {
      return new SimpleSequenceFactory().createSequence(
        createDNA(dna),
        "", name, new SimpleAnnotation()
      );
    //} catch (BioException se) {
    //  throw new BioError("Something has gone badly wrong with DNA", se);
    //}
  }


    /** Get a new dna as a GappedSequence */
    public static GappedSequence createGappedDNASequence(String dna, String name) throws IllegalSymbolException{
        String dna1 = dna.replaceAll("-", "");
        Sequence dnaSeq = createDNASequence(dna1, name);
        GappedSequence dnaSeq1 = new SimpleGappedSequence(dnaSeq);
        int pos = dna.indexOf('-', 0);
        while(pos!=-1){
            dnaSeq1.addGapInView(pos+1);
            pos = dna.indexOf('-', pos+1);
        }
        return dnaSeq1;
    }


  /**
   * Return an integer index for a symbol - compatible with
   * <code>forIndex</code>.
   *
   * <p>
   * The index for a symbol is stable accross virtual machines &
   * invocations.
   * </p>
   *
   * @param sym  the Symbol to index
   * @return the index for that symbol
   *
   * @throws IllegalSymbolException if sym is not a member of the DNA
   * alphabet
   */
  public static int index(Symbol sym) throws IllegalSymbolException {
    if(sym == a) {
      return 0;
    } else if(sym == g) {
      return 1;
    } else if(sym == c) {
      return 2;
    } else if(sym == t) {
      return 3;
    }
    getDNA().validate(sym);
    throw new IllegalSymbolException("Really confused. Can't find index for " +
                                      sym.getName());
  }

  /**
   * Return the symbol for an index - compatible with <code>index</code>.
   *
   * <p>
   * The index for a symbol is stable accross virtual machines &
   * invocations.
   * </p>
   *
   * @param index  the index to look up
   * @return       the symbol at that index
   *
   * @throws IndexOutOfBoundsException if index is not between 0 and 3
   */
  static public Symbol forIndex(int index)
  throws IndexOutOfBoundsException {
    if(index == 0)
      return a;
    else if(index == 1)
      return g;
    else if(index == 2)
      return c;
    else if(index == 3)
      return t;
    else throw new IndexOutOfBoundsException("No symbol for index " + index);
  }

  /**
   * Complement the symbol.
   *
   * @param sym  the symbol to complement
   * @return a Symbol that is the complement of sym
   * @throws IllegalSymbolException if sym is not a member of the DNA alphabet
   */
  static public Symbol complement(Symbol sym)
  throws IllegalSymbolException {
    if(sym == a) {
      return t;
    } else if(sym == g) {
      return c;
    } else if(sym == c) {
      return g;
    } else if(sym == t) {
      return a;
    }
    Symbol s = symbolToComplement.get(sym);
    if(s != null) {
      return s;
    } else {
      getDNA().validate(sym);
      throw new BioError(
        "Really confused. Can't find symbol " +
        sym.getName()
      );
    }
  }

  /**
   * Retrieve the symbol for a symbol.
   *
   * @param token  the char to look up
   * @return  the symbol for that char
   * @throws IllegalSymbolException if the char is not a valid IUB dna code
   */
  static public Symbol forSymbol(char token)
  throws IllegalSymbolException {
    String t = String.valueOf(token);
    SymbolTokenization toke;
    try{
      toke = getDNA().getTokenization("token");
    }catch(BioException e){
      throw new BioError("Cannot find the 'token' Tokenization for DNA!?", e);
    }
    return toke.parseToken(t);
  }

  /**
   * Retrieve a complement view of list.
   *
   * @param list  the SymbolList to complement
   * @return a SymbolList that is the complement
   * @throws IllegalAlphabetException if list is not a complementable alphabet
   */
  public static SymbolList complement(SymbolList list)
  throws IllegalAlphabetException {
    return SymbolListViews.translate(list, complementTable());
  }

  /**
   * Retrieve a reverse-complement view of list.
   *
   * @param list  the SymbolList to complement
   * @return a SymbolList that is the complement
   * @throws IllegalAlphabetException if list is not a complementable alphabet
   */
  public static SymbolList reverseComplement(SymbolList list)
  throws IllegalAlphabetException {
    return SymbolListViews.translate(SymbolListViews.reverse(list), complementTable());
  }

  /**
   * Returns a SymbolList that is reverse complemented if the strand is
   * negative, and the origninal one if it is not.
   *
   * @param list  the SymbolList to view
   * @param strand the Strand to use
   * @return the apropreate view of the SymbolList
   * @throws IllegalAlphabetException if list is not a complementable alphabet
   */
  public static SymbolList flip(SymbolList list, StrandedFeature.Strand strand)
  throws IllegalAlphabetException {
    if(strand == StrandedFeature.NEGATIVE) {
      return reverseComplement(list);
    } else {
      return list;
    }
  }

  /**
   * Get a translation table for complementing DNA symbols.
   *
   * @since 1.1
   */

  public static ReversibleTranslationTable complementTable() {
    return complementTable;
  }

    /**
     * Get a single-character token for a DNA symbol
     *
     * @throws IllegalSymbolException if <code>sym</code> is not a member of the DNA alphabet
     */

    public static char dnaToken(Symbol sym)
        throws IllegalSymbolException
    {
        return dnaTokens.tokenizeSymbol(sym).charAt(0);
    }

  /**
   * Sneaky class for complementing DNA bases.
   */

  private static class DNAComplementTranslationTable
  implements ReversibleTranslationTable {
    public Symbol translate(Symbol s)
          throws IllegalSymbolException {
            return DNATools.complement(s);
          }

    public Symbol untranslate(Symbol s)
          throws IllegalSymbolException {
            return DNATools.complement(s);
          }

          public Alphabet getSourceAlphabet() {
            return DNATools.getDNA();
          }

          public Alphabet getTargetAlphabet() {
            return DNATools.getDNA();
          }
  }

  /**
   * return a SimpleDistribution of specified GC content.
   * @param fractionGC (G+C) content as a fraction.
   */
  public static Distribution getDNADistribution(double fractionGC)
  {
    try {
        Distribution dist = new SimpleDistribution(DNATools.getDNA());
        double gc = 0.5 * fractionGC;
        double at = 0.5 * (1.0 - fractionGC);
        dist.setWeight(DNATools.a(), at);
        dist.setWeight(DNATools.t(), at);
        dist.setWeight(DNATools.c(), gc);
        dist.setWeight(DNATools.g(), gc);

        return dist;
    }
        // these exceptions are just plain impossible!!!
    catch (IllegalSymbolException ise) { return null; }
    catch (ChangeVetoException cve) { return null; }
  }

  /**
   * return a (DNA x DNA) cross-product Distribution with specified
   * DNA contents in each component Alphabet.
   * @param fractionGC0 (G+C) content of first sequence as a fraction.
   * @param fractionGC1 (G+C) content of second sequence as a fraction.
   */
  public static Distribution getDNAxDNADistribution(
    double fractionGC0,
    double fractionGC1
    )
  {
    return new PairDistribution(getDNADistribution(fractionGC0), getDNADistribution(fractionGC1));
  }
  
  /**
   * Converts a <code>SymbolList</code> from the DNA <code>Alphabet</code> to the
   * RNA <code>Alphabet</code>.
   * @param syms the <code>SymbolList</code> to convert to RNA
   * @return a view on <code>syms</code> where <code>Symbols</code> have been converted to RNA.
   * Most significantly t's are now u's. The 5' to 3' order of the Symbols is conserved.
   * @since 1.4
   * @throws IllegalAlphabetException if <code>syms</code> is not DNA.
   */
   public static SymbolList toRNA(SymbolList syms)throws IllegalAlphabetException{
     return SymbolListViews.translate(syms, transcriptionTable);
   }
   
   /**
    * Transcribes DNA to RNA. The method more closely represents the biological reality
    * than <code>toRNA(SymbolList syms)</code> does. The presented DNA <code>SymbolList</code> 
    * is assumed to be the template strand in the 5' to 3' orientation. The resulting
    * RNA is transcribed from this template effectively a reverse complement in the RNA alphabet.
    * The method is equivalent to calling <code>reverseComplement()</code> and <code>toRNA()</code> in sequence.
    * <p>If you are dealing with cDNA sequences that you want converted to RNA you would be 
    * better off calling <code>toRNA(SymbolList syms)</code>
    * @param syms the <code>SymbolList</code> to convert to RNA
    * @return a view on <code>syms</code> where <code>Symbols</code> have been converted to RNA.
    * @since 1.4
    * @throws IllegalAlphabetException if <code>syms</code> is not DNA.
    */
   public static SymbolList transcribeToRNA(SymbolList syms) throws IllegalAlphabetException{
     syms = reverseComplement(syms);
     return toRNA(syms);
   }
   
   /**
    * Convenience method that directly converts a DNA sequence to RNA then to
    * protein. The translated protein is from the +1 reading frame of the
    * <code>SymbolList</code>. The whole <code>SymbolList</code> is translated
    * although up to 2 DNA residues may be truncated if full codons cannot be 
    * formed.
    * @param syms the sequence to be translated.
    * @return the translated protein sequence.
    * @throws org.biojava.bio.symbol.IllegalAlphabetException if <code>syms</code>
    * is not from the DNA alphabet.
    * @since 1.5.1
    */
   public static SymbolList toProtein(final SymbolList syms) 
           throws IllegalAlphabetException{
       SymbolList symz = new SimpleSymbolList(syms);
       symz = toRNA(symz);
       //truncate to a length divisible by three.
       symz = symz.subList(1, symz.length() - (symz.length() %3));
       return RNATools.translate(symz);
   }
   
   /**
    * Convenience method to translate a region of a DNA sequence directly into
    * protein. While the start and end can be specified if the length of the 
    * specified region is not evenly divisible by three then the translated 
    * region will be truncated until a full terminal codon can be formed.
    * @param syms the DNA sequence to be translated.
    * @param start the location to begin translation.
    * @param end the end of the translated region.
    * @return the translated protein sequence.
    * @throws org.biojava.bio.symbol.IllegalAlphabetException if <code>syms
    * </code> is not from the DNA alphabet.
    * @since 1.5.1
    */
   public static SymbolList toProtein(final SymbolList syms, int start, int end)
           throws IllegalAlphabetException{
       SymbolList symz = new SimpleSymbolList(syms);
       symz = symz.subList(start, end);
       //return toProtein(syms);
       //changed to fix bug 2521
       return toProtein(symz);
   }
}






