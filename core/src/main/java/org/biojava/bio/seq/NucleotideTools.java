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
import org.biojava.bio.seq.impl.SimpleSequenceFactory;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.ReversibleTranslationTable;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;



/**

 * Useful functionality for processing nucleotide sequences.

 *

 * @author Matthew Pocock

 * @author Keith James (docs)

 */

public final class NucleotideTools {

  private static final ReversibleTranslationTable complementTable;

  static private final FiniteAlphabet nucleotide;

    private static final SymbolTokenization nucleotideTokens;



  static private final AtomicSymbol a;

  static private final AtomicSymbol g;

  static private final AtomicSymbol c;

  static private final AtomicSymbol t;

  static private final AtomicSymbol u;

  static private final Symbol r;

  static private final Symbol y;

  static private final Symbol m;

  static private final Symbol k;

  static private final Symbol s;

  static private final Symbol w;

  static private final Symbol b;

  static private final Symbol d;

  static private final Symbol h;

  static private final Symbol v;

  static private final Symbol n;





  static private Map symbolToComplement;



  static {

    try {

      nucleotide = (FiniteAlphabet) AlphabetManager.alphabetForName("NUCLEOTIDE");

      nucleotideTokens = nucleotide.getTokenization("token");

      SymbolList syms = new SimpleSymbolList(nucleotideTokens, "agcturymkswbdhvn");

      a = (AtomicSymbol) syms.symbolAt(1);

      g = (AtomicSymbol) syms.symbolAt(2);

      c = (AtomicSymbol) syms.symbolAt(3);

      t = (AtomicSymbol) syms.symbolAt(4);

      u = (AtomicSymbol) syms.symbolAt(5);

      r = syms.symbolAt(6);

      y = syms.symbolAt(7);

      m = syms.symbolAt(8);

      k = syms.symbolAt(9);

      s = syms.symbolAt(10);

      w = syms.symbolAt(11);

      b = syms.symbolAt(12);

      d = syms.symbolAt(13);

      h = syms.symbolAt(14);

      v = syms.symbolAt(15);

      n = syms.symbolAt(16);



      symbolToComplement = new HashMap();



      // add the gap symbol

      Symbol gap = nucleotide.getGapSymbol();

      symbolToComplement.put(gap, gap);



      // add all other ambiguity symbols

      for(Iterator i = AlphabetManager.getAllSymbols(nucleotide).iterator(); i.hasNext();) {

          Symbol as = (Symbol) i.next();

          FiniteAlphabet matches = (FiniteAlphabet) as.getMatches();

          if (matches.size() > 1) {   // We've hit an ambiguous symbol.

              Set l = new HashSet();

              for(Iterator j = matches.iterator(); j.hasNext(); ) {

                  l.add(complement((Symbol) j.next()));

              }

              symbolToComplement.put(as, nucleotide.getAmbiguity(l));

          }

      }





      complementTable = new NucleotideComplementTranslationTable();

    } catch (Throwable t) {

      throw new BioError("Unable to initialize NucleotideTools",t);

    }

  }



  public static AtomicSymbol a() { return a; }

  public static AtomicSymbol g() { return g; }

  public static AtomicSymbol c() { return c; }

  public static AtomicSymbol t() { return t; }

  public static AtomicSymbol u() { return u; }

  public static Symbol r() { return r; }

  public static Symbol y() { return y; }

  public static Symbol m() { return m; }

  public static Symbol k() { return k; }

  public static Symbol s() { return s; }

  public static Symbol w() { return w; }

  public static Symbol b() { return b; }

  public static Symbol d() { return d; }

  public static Symbol h() { return h; }

  public static Symbol v() { return v; }

  public static Symbol n() { return n; }


  private NucleotideTools() {
  }

  /**

   * Return the Nucleotide alphabet.

   *

   * @return a flyweight version of the Nucleotide alphabet

   */

  public static FiniteAlphabet getNucleotide() {

    return nucleotide;

  }



  /**

   * Return a new Nucleotide <span class="type">SymbolList</span> for

   * <span class="arg">nucleotide</span>.

   *

   * @param nucleotide a <span class="type">String</span> to parse into Nucleotide

   * @return a <span class="type">SymbolList</span> created form

   *         <span class="arg">nucleotide</span>

   * @throws IllegalSymbolException if <span class="arg">nucleotide</span> contains

   *         any non-Nucleotide characters

   */

  public static SymbolList createNucleotide(String nucleotide)

  throws IllegalSymbolException {

    try {

      SymbolTokenization p = getNucleotide().getTokenization("token");

      return new SimpleSymbolList(p, nucleotide);

    } catch (BioException se) {

      throw new BioError("Something has gone badly wrong with Nucleotide",se);

    }

  }



  /**

   * Return a new Nucleotide <span class="type">Sequence</span> for

   * <span class="arg">nucleotide</span>.

   *

   * @param nucleotide a <span class="type">String</span> to parse into Nucleotide

   * @param name a <span class="type">String</span> to use as the name

   * @return a <span class="type">Sequence</span> created form

   *         <span class="arg">nucleotide</span>

   * @throws IllegalSymbolException if <span class="arg">nucleotide</span> contains

   *         any non-Nucleotide characters

   */

  public static Sequence createNucleotideSequence(String nucleotide, String name)

  throws IllegalSymbolException {

    try {

      return new SimpleSequenceFactory().createSequence(

        createNucleotide(nucleotide),

        "", name, new SimpleAnnotation()

      );

    } catch (BioException se) {

      throw new BioError("Something has gone badly wrong with Nucleotide",se);

    }

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

   * @throws IllegalSymbolException if sym is not a member of the Nucleotide

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

    } else if(sym == u) {

      return 4;

    }

    getNucleotide().validate(sym);

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

    else if(index == 4)

      return u;

    else throw new IndexOutOfBoundsException("No symbol for index " + index);

  }



  /**

   * Complement the symbol.

   *

   * @param sym  the symbol to complement

   * @return a Symbol that is the complement of sym

   * @throws IllegalSymbolException if sym is not a member of the Nucleotide alphabet

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

    } else if(sym == u) {

      return a;

    }

    Symbol s = (Symbol) symbolToComplement.get(sym);

    if(s != null) {

      return s;

    } else {

      getNucleotide().validate(sym);

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

   * @throws IllegalSymbolException if the char does not belong to {a, g, c, t, u}

   */

  static public Symbol forSymbol(char token)

  throws IllegalSymbolException {

    if(token == 'a') {

      return a;

    } else if(token == 'g') {

      return g;

    } else if(token == 'c') {

      return c;

    } else if(token == 't') {

      return t;

    } else if(token == 'u') {

      return u;

    }

    throw new IllegalSymbolException("Unable to find symbol for token " + token);

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

   * Get a translation table for complementing Nucleotide symbols.

   *

   * @since 1.1

   */



  public static ReversibleTranslationTable complementTable() {

    return complementTable;

  }



    /**

     * Get a single-character token for a Nucleotide symbol

     *

     * @throws IllegalSymbolException if <code>sym</code> is not a member of the Nucleotide alphabet

     */



    public static char nucleotideToken(Symbol sym)

        throws IllegalSymbolException

    {

        return nucleotideTokens.tokenizeSymbol(sym).charAt(0);

    }



  /**

   * Sneaky class for complementing Nucleotide bases.

   */



  private static class NucleotideComplementTranslationTable

  implements ReversibleTranslationTable {

    public Symbol translate(Symbol s)

          throws IllegalSymbolException {

            return NucleotideTools.complement(s);

          }



    public Symbol untranslate(Symbol s)

          throws IllegalSymbolException	{

            return NucleotideTools.complement(s);

          }



          public Alphabet getSourceAlphabet() {

            return NucleotideTools.getNucleotide();

          }



          public Alphabet getTargetAlphabet() {

            return NucleotideTools.getNucleotide();

          }

  }

}



