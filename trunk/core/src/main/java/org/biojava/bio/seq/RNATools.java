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

import java.io.InputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.impl.SimpleSequenceFactory;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AbstractReversibleTranslationTable;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.ManyToOneTranslationTable;
import org.biojava.bio.symbol.ReversibleTranslationTable;
import org.biojava.bio.symbol.SimpleGeneticCodeTable;
import org.biojava.bio.symbol.SimpleReversibleTranslationTable;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.utils.ClassTools;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

/**
 * Useful functionality for processing DNA and RNA sequences.
 *
 * @author Matthew Pocock
 * @author Keith James (docs)
 * @author Thomas Down
 * @author Greg Cox
 * @author Mark Schreiber
 * @author David Huen (refactoring)
 * @author gwaldon (update genetic code translation tables)
 */
public final class RNATools {
  private static final ReversibleTranslationTable complementTable;
  private static final SimpleReversibleTranslationTable transcriptionTable;
  static private final FiniteAlphabet rna;
  static private final Map geneticCodes;

  static private final AtomicSymbol a;
  static private final AtomicSymbol g;
  static private final AtomicSymbol c;
  static private final AtomicSymbol u;
  static private final Symbol n;

  static private Map symbolToComplement;

  static {
    try {
      rna = (FiniteAlphabet) AlphabetManager.alphabetForName("RNA");

      SymbolList syms = new SimpleSymbolList(rna.getTokenization("token"), "agcun");
      a = (AtomicSymbol) syms.symbolAt(1);
      g = (AtomicSymbol) syms.symbolAt(2);
      c = (AtomicSymbol) syms.symbolAt(3);
      u = (AtomicSymbol) syms.symbolAt(4);
      n = syms.symbolAt(5);

      symbolToComplement = new HashMap();

      // add the gap symbol
      Symbol gap = rna.getGapSymbol();
      symbolToComplement.put(gap, gap);

      // add all other ambiguity symbols
      for(Iterator i = AlphabetManager.getAllSymbols(rna).iterator(); i.hasNext();) {
          Symbol as = (Symbol) i.next();
          FiniteAlphabet matches = (FiniteAlphabet) as.getMatches();
          if (matches.size() > 1) {   // We've hit an ambiguous symbol.
              Set l = new HashSet();
              for(Iterator j = matches.iterator(); j.hasNext(); ) {
                  l.add(complement((Symbol) j.next()));
              }
              symbolToComplement.put(as, rna.getAmbiguity(l));
          }
      }
      complementTable = new RNAComplementTranslationTable();

      transcriptionTable = new SimpleReversibleTranslationTable(DNATools.getDNA(), rna);
      transcriptionTable.setTranslation(DNATools.a(), a);
      transcriptionTable.setTranslation(DNATools.c(), c);
      transcriptionTable.setTranslation(DNATools.g(), g);
      transcriptionTable.setTranslation(DNATools.t(), u);

      geneticCodes = new HashMap();
      loadGeneticCodes();
    } catch (Throwable t) {
      throw new BioError("Unable to initialize RNATools", t);
    }
  }

  public static AtomicSymbol a() { return a; }
  public static AtomicSymbol g() { return g; }
  public static AtomicSymbol c() { return c; }
  public static AtomicSymbol u() { return u; }
  public static Symbol n() { return n; }

  private RNATools() {
  }
  
  /**
   * Return the RNA alphabet.
   *
   * @return a flyweight version of the RNA alphabet
   */
  public static FiniteAlphabet getRNA() {
    return rna;
  }

  /**
   * Gets the (RNA x RNA x RNA) Alphabet
   * @return a flyweight version of the (RNA x RNA x RNA) alphabet
   */
  public static FiniteAlphabet getCodonAlphabet(){
    return (FiniteAlphabet)AlphabetManager.generateCrossProductAlphaFromName("(RNA x RNA x RNA)");
  }

  /**
   * Return a new RNA <span class="type">SymbolList</span> for
   * <span class="arg">rna</span>.
   *
   * @param rna a <span class="type">String</span> to parse into RNA
   * @return a <span class="type">SymbolList</span> created form
   *         <span class="arg">rna</span>
   * @throws IllegalSymbolException if  <span class="arg">rna</span> contains
   *         any non-RNA characters
   */
  public static SymbolList createRNA(String rna)
  throws IllegalSymbolException {
    SymbolTokenization p = null;
    try {
      p = getRNA().getTokenization("token");
    } catch (BioException e) {
      throw new BioError("Something has gone badly wrong with RNA", e);
    }
    return new SimpleSymbolList(p, rna);

  }

  /**
   * Return a new RNA <span class="type">Sequence</span> for
   * <span class="arg">rna</span>.
   *
   * @param rna a <span class="type">String</span> to parse into RNA
   * @param name a <span class="type">String</span> to use as the name
   * @return a <span class="type">Sequence</span> created form
   *         <span class="arg">dna</span>
   * @throws IllegalSymbolException if <span class="arg">rna</span> contains
   *         any non-DNA characters
   */
  public static Sequence createRNASequence(String rna, String name)
  throws IllegalSymbolException {
    try {
      return new SimpleSequenceFactory().createSequence(
        createRNA(rna),
        "", name, new SimpleAnnotation()
      );
    } catch (BioException se) {
      throw new BioError("Something has gone badly wrong with RNA", se);
    }
  }

  /**
   * Return an integer index for a symbol - compatible with forIndex.
   * <p>
   * The index for a symbol is stable across virtual machines & invocations.
   *
   * @param sym  the Symbol to index
   * @return     the index for that symbol
   * @throws IllegalSymbolException if sym is not a member of the DNA alphabet
   */
  public static int index(Symbol sym) throws IllegalSymbolException {
    if(sym == a) {
      return 0;
    } else if(sym == g) {
      return 1;
    } else if(sym == c) {
      return 2;
    } else if(sym == u) {
      return 3;
    }
    getRNA().validate(sym);
    throw new IllegalSymbolException("Really confused. Can't find index for " +
                                      sym.getName());
  }

  /**
   * Return the symbol for an index - compatible with index.
   * <p>
   * The index for a symbol is stable accross virtual machines & invocations.
   *
   * @param index  the index to look up
   * @return       the symbol at that index
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
      return u;
    else throw new IndexOutOfBoundsException("No symbol for index " + index);
  }

  /**
   * Complement the symbol.
   *
   * @param sym  the symbol to complement
   * @return a Symbol that is the complement of sym
   * @throws IllegalSymbolException if sym is not a member of the RNA alphabet
   */
  static public Symbol complement(Symbol sym)
  throws IllegalSymbolException {
    if(sym == a) {
      return u;
    } else if(sym == g) {
      return c;
    } else if(sym == c) {
      return g;
    } else if(sym == u) {
      return a;
    }
    Symbol s = (Symbol) symbolToComplement.get(sym);
    if(s != null) {
      return s;
    } else {
      getRNA().validate(sym);
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
   * @throws IllegalSymbolException if the char is not a valid IUB code.
   */
  static public Symbol forSymbol(char token)
  throws IllegalSymbolException {
    String t = String.valueOf(token);
    SymbolTokenization toke;

    try{
      toke = getRNA().getTokenization("token");
    }catch(BioException e){
      throw new BioError("Cannot find the 'token' Tokenization for RNA!?", e);
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
   * Transcribe DNA into RNA.
   * @deprecated The naming of this method is confusing and inconsistent use either DNATools.toRNA(SymbolList list) or
   * DNATools.transcribeToRNA(SymbolList list) depending on the desired behaivour.
   * @param list the SymbolList to transcribe
   * @return a SymbolList that is the transcribed view
   * @throws IllegalAlphabetException if the list is not DNA
   */
   public static SymbolList transcribe(SymbolList list)
   throws IllegalAlphabetException {
     return SymbolListViews.translate(list, transcriptionTable());
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
   * Get a translation table for converting DNA to RNA.
   *
   * @since 1.1
   */
  public static ReversibleTranslationTable transcriptionTable() {
    return transcriptionTable;
  }

  /**
   * Retrieve a TranslationTable by name. The valid names are:
   *
   * <ul>
   * <li>"UNIVERSAL"</li>
   * <li>"VERTEBRATE_MITOCHONDRIAL"</li>
   * <li>"YEAST_MITOCHONDRIAL"</li>
   * <li>"MOLD_MITOCHONDRIAL"</li>
   * <li>"INVERTEBRATE_MITOCHONDRIAL"</li>
   * <li>"CILIATE_NUCLEAR"</li>
   * <li>"ECHINODERM_MITOCHONDRIAL"</li>
   * <li>"EUPLOTID_NUCLEAR"</li>
   * <li>"BACTERIAL"</li>
   * <li>"ALTERNATIVE_YEAST_NUCLEAR"</li>
   * <li>"ASCIDIAN_MITOCHONDRIAL"</li>
   * <li>"FLATWORM_MITOCHONDRIAL"</li>
   * <li>"BLEPHARISMA_MACRONUCLEAR"</li>
   * <li>"CHLOROPHYCEAN_MITOCHONDRIAL"</li>
   * <li>"TREMATODE_MITOCHONDRIAL"</li>
   * <li>"SCENEDESMUS_MITOCHONDRIAL"</li>
   * </ul>
   *
   * There are public static final fields in the TranslationTable
   * interface which contain these values. One of these should be used
   * as the argument for this method.
   * <p>
   * You can now get the reverse translation of the residue back to its
   * (usually several) codons too.
   *
   * @since 1.1
   */
  public static ManyToOneTranslationTable getGeneticCode(String name) {
    return (ManyToOneTranslationTable) geneticCodes.get(name);
  }
  
  /**
   * Retrieve a TranslationTable by number.
   * These numbers correspond to the transl_table qualifier in the
   * DDBJ/EMBL/GenBank Feature Table (Version 6.5  Apr 2006): transl_table
   * defines the genetic code table used if other than the universal 
   * genetic code table. Tables are described in appendix V,
   * section 7.5.5:
   *
   * <ul>
   * <li>" 1 - UNIVERSAL"</li>
   * <li>" 2 - VERTEBRATE_MITOCHONDRIAL"</li>
   * <li>" 3 - YEAST_MITOCHONDRIAL"</li>
   * <li>" 4 - MOLD_MITOCHONDRIAL"</li>
   * <li>" 5 - INVERTEBRATE_MITOCHONDRIAL"</li>
   * <li>" 6 - CILIATE_NUCLEAR"</li>
   * <li>" 9 - ECHINODERM_MITOCHONDRIAL"</li>
   * <li>"10 - EUPLOTID_NUCLEAR"</li>
   * <li>"11 - BACTERIAL"</li>
   * <li>"12 - ALTERNATIVE_YEAST_NUCLEAR"</li>
   * <li>"13 - ASCIDIAN_MITOCHONDRIAL"</li>
   * <li>"14 - FLATWORM_MITOCHONDRIAL"</li>
   * <li>"15 - BLEPHARISMA_MACRONUCLEAR"</li>
   * <li>"16 - 2CHLOROPHYCEAN_MITOCHONDRIAL"</li>
   * <li>"21 - TREMATODE_MITOCHONDRIAL"</li>
   * <li>"23 - SCENEDESMUS_MITOCHONDRIAL"</li>
   * </ul>
   *
   * @throws IllegalArgumentException if there is no table with that number.
   * @since 1.5
   */
  public static ManyToOneTranslationTable getGeneticCode(int table_num) {
      Set tables = getGeneticCodeNames();
      Iterator it = tables.iterator();
      while(it.hasNext()) {
          String tableName = (String) it.next();
          SimpleGeneticCodeTable table = (SimpleGeneticCodeTable) geneticCodes.get(tableName);
          if(table.getTableNumber()==table_num)
              return table;
      }
      throw new IllegalArgumentException("There is no genetic code table at that number");
  }
  
  /**
   * Retrieve a Set containing the name of each genetic code.
   *
   * @since 1.1
   */
  public static Set getGeneticCodeNames() {
    return geneticCodes.keySet();
  }

  /**
   * Translate RNA into protein (with termination symbols).  For
   * compatibility with BioJava 1.1, this will also handle sequences
   * which are already expressed in the (RNA x RNA x RNA) (codon)
   * alphabet.
   *
   * @since 1.1
   */
  public static SymbolList translate(SymbolList syms)
    throws IllegalAlphabetException
  {
      if (syms.getAlphabet() == getRNA()) {
          syms = SymbolListViews.windowedSymbolList(syms, 3);
      }
      return SymbolListViews.translate(syms, getGeneticCode("UNIVERSAL"));
  }

  private static void loadGeneticCodes() {
    try {
      InputStream tablesStream = ClassTools.getClassLoader(RNATools.class).getResourceAsStream(
        "org/biojava/bio/seq/TranslationTables.xml"
      );
      if(tablesStream == null ) {
        throw new BioError("Couldn't locate TranslationTables.xml.");
      }

      InputSource is = new InputSource(tablesStream);
      DocumentBuilder parser = DocumentBuilderFactory.newInstance().newDocumentBuilder();
      Document doc = parser.parse(is);

      NodeList children = doc.getDocumentElement().getChildNodes();
      for(int i = 0; i < children.getLength(); i++) {
        Node cnode = children.item(i);
        if(! (cnode instanceof Element)) {
          continue;
        }

        Element child = (Element) cnode;
        String name = child.getNodeName();
        if(name.equals("table")) {
          String tableName = child.getAttribute("name");
          String source = child.getAttribute("source");
          String target = child.getAttribute("target");
          FiniteAlphabet sourceA =
            (FiniteAlphabet) AlphabetManager.alphabetForName(source);
          FiniteAlphabet targetA =
            (FiniteAlphabet) AlphabetManager.alphabetForName(target);
          SymbolTokenization targetP = targetA.getTokenization("name");
          SimpleGeneticCodeTable table = new SimpleGeneticCodeTable (
            sourceA,
            targetA
          );

          NodeList translates = child.getChildNodes();
          for(int j = 0; j < translates.getLength(); j++) {
            Node tn = translates.item(j);
            if(tn instanceof Element) {
              Element te = (Element) tn;
              if(te.getTagName().equals("transl_table")) {
                  int num = Integer.valueOf(te.getAttribute("value")).intValue();
                  String description = te.getAttribute("description");
                  table.setTableNumber(num);
                  table.setDescription(description);
                  continue;
              }
              String from = te.getAttribute("from");
              String to = te.getAttribute("to");

              //
              // Not the most elegant solution, but I wanted this working
              // quickly for 1.1.  It's been broken for ages.
              //     -td 26/i/20001
              //

              SymbolList fromSymbols = RNATools.createRNA(from);
              if (fromSymbols.length() != 3) {
                  throw new BioError("`" + from + "' is not a valid codon");
              }

              // AtomicSymbol fromS = (AtomicSymbol) sourceP.parseToken(from);
              AtomicSymbol fromS = (AtomicSymbol) sourceA.getSymbol(fromSymbols.toList());
              AtomicSymbol toS   = (AtomicSymbol) targetP.parseToken(to);
              table.setTranslation(fromS, toS);
            }
          }

          geneticCodes.put(tableName, table);
        }
      }
    } catch (Exception e) {
      throw new BioError("Couldn't parse TranslationTables.xml", e);
    }
  }

  /**
   * Sneaky class for complementing RNA bases.
   */

  private static class RNAComplementTranslationTable
  extends AbstractReversibleTranslationTable {
    public Symbol doTranslate(Symbol s)
          throws IllegalSymbolException {
            return (Symbol) RNATools.complement(s);
          }

          public Symbol doUntranslate(Symbol s)
          throws IllegalSymbolException {
            return (Symbol) RNATools.complement(s);
    }

          public Alphabet getSourceAlphabet() {
            return RNATools.getRNA();
          }

          public Alphabet getTargetAlphabet() {
            return RNATools.getRNA();
          }
  }
}

