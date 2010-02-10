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


package org.biojava.bio.symbol;

import java.io.IOException;
import java.io.InputStream;
import java.io.InvalidObjectException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.WeakHashMap;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.io.AlternateTokenization;
import org.biojava.bio.seq.io.CharacterTokenization;
import org.biojava.bio.seq.io.NameTokenization;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.StreamParser;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ClassTools;
import org.biojava.utils.Unchangeable;
import org.biojava.utils.cache.WeakValueHashMap;
import org.biojava.utils.lsid.Identifiable;
import org.biojava.utils.lsid.LifeScienceIdentifier;
import org.biojava.utils.lsid.LifeScienceIdentifierParseException;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.SAX2StAXAdaptor;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;


/**
 * Utility methods for working with Alphabets.  Also acts as a registry for
 * well-known alphabets.
 *
 * <p>
 * The alphabet interfaces themselves don't give you a lot of help in actually
 * getting an alphabet instance. This is where the AlphabetManager comes in
 * handy. It helps out in serialization, generating derived alphabets and
 * building CrossProductAlphabet instances. It also contains limited support for
 * parsing complex alphabet names back into the alphabets.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Mark Schreiber
 * @author George Waldon (alternate tokenization)
 */

public final class AlphabetManager {
  static private Map nameToAlphabet;
  //static private Map nameToSymbol;
  static private Map lsidToSymbol;
  static private Map crossProductAlphabets;
  static private Map ambiguitySymbols;
  static private GapSymbol gapSymbol;
  static private Map gapBySize;
  static private Map alphabetToIndex = new WeakHashMap();
  static private Map symListToSymbol;
    
  /**
   * <p>
   * Initialize the static AlphabetManager resources.
   * </p>
   *
   * <p>
   * This parses the resource
   * <code>org/biojava/bio/seq/tools/AlphabetManager.xml</code>
   * and builds a basic set of alphabets.
   * </p>
   */
  static {
    nameToAlphabet = new HashMap();
    //nameToSymbol = new HashMap();
    lsidToSymbol = new HashMap();
    ambiguitySymbols = new HashMap();

    gapSymbol = new GapSymbol();
    gapBySize = new HashMap();
    gapBySize.put(new SizeQueen(new ArrayList()), gapSymbol);

    nameToAlphabet.put("INTEGER", IntegerAlphabet.getInstance());
    nameToAlphabet.put("DOUBLE", DoubleAlphabet.getInstance());

    symListToSymbol = new WeakValueHashMap();

    try {
      SizeQueen sq = new SizeQueen(Arrays.asList(
                new Alphabet[] { DoubleAlphabet.getInstance() }));  
      gapBySize.put(sq, 
                    new WellKnownGapSymbol(
                         Arrays.asList(new Symbol[] { gapSymbol}), sq));
    } catch (IllegalSymbolException ise) {
      throw new BioError(

        "Assertion Failure: Should be able to make gap basis", ise
      );
    }

    ambiguitySymbols.put(new HashSet(), gapSymbol);
    try {
      InputStream alphabetStream = ClassTools.getClassLoader(AlphabetManager.class).getResourceAsStream(
        "org/biojava/bio/symbol/AlphabetManager.xml"
      );
      if (alphabetStream == null) {
          throw new BioError("Couldn't locate AlphabetManager.xml.  This probably means that your biojava.jar file is corrupt or incorrectly built.");
      }
      InputSource is = new InputSource(alphabetStream);
      loadAlphabets(is);
    } catch (Exception t) {
      throw new BioError( "Unable to initialize AlphabetManager", t);
    }
  }

    /**
   * Singleton instance.
   */
  static private AlphabetManager am;

  /**
   * Retrieve the singleton instance.
   *
   * @return the AlphabetManager instance
   * @deprecated all AlphabetManager methods have become static
   */
  static public AlphabetManager instance() {
    if(am == null)
      am = new AlphabetManager();
    return am;
  }

  
    /**
     * Return the ambiguity symbol which matches all symbols in
     * a given alphabet.
     * @since 1.2
     * @param alpha The alphabet
     * @return the ambiguity symbol
     */

    public static Symbol getAllAmbiguitySymbol(FiniteAlphabet alpha) {
        Set allSymbols = new HashSet();
        for (Iterator i = alpha.iterator(); i.hasNext(); ) {
            allSymbols.add(i.next());
        }
        try {
            return alpha.getAmbiguity(allSymbols);
        } catch (IllegalSymbolException ex) {
            throw new BioError( "Assertion failure: coudn't recover all-ambiguity symbol", ex);
        }
    }

    /**
     * Return a set containing all possible symbols which can be
     * considered members of a given alphabet, including ambiguous
     * symbols.  Warning, this method can return large sets!
     * @since 1.2
     * @param alpha The alphabet
     * @return The set of symbols that are members of <code>alpha</code>
     */

    public static Set getAllSymbols(FiniteAlphabet alpha) {
        Set allSymbols = new HashSet();
        List orderedAlpha = new ArrayList(alpha.size());
        for (Iterator i = alpha.iterator(); i.hasNext(); ) {
            orderedAlpha.add(i.next());
        }

        int atomicSyms = alpha.size();
        int totalSyms = 1 << atomicSyms;

        for (int cnt = 0; cnt < totalSyms; ++cnt) {
            Set matchSet = new HashSet();
            for (int atom = 0; atom < atomicSyms; ++atom) {
                if ((cnt & (1 << atom)) != 0) {
                    matchSet.add(orderedAlpha.get(atom));
                }
            }

            try {
                allSymbols.add(alpha.getAmbiguity(matchSet));
            } catch (IllegalSymbolException ex) {
                throw new BioError( "Assertion failed: couldn't get ambiguity symbol", ex);
            }
        }

        return allSymbols;
    }



  /**
   * Retrieve the alphabet for a specific name.
   *
   * @param name the name of the alphabet
   * @return the alphabet object
   * @throws NoSuchElementException if there is no alphabet by that name
   */
  static public Alphabet alphabetForName(String name)
  throws NoSuchElementException{
    Alphabet alpha = (Alphabet) nameToAlphabet.get(name);
    if(alpha == null) {
      if(name.startsWith("(") && name.endsWith(")")) {
        alpha = generateCrossProductAlphaFromName(name);
      } else {
        throw new NoSuchElementException(
          "No alphabet for name " + name + " could be found"
        );
      }
    }
    return alpha;
  }
   /**
   * Retrieve the symbol represented a String object
   * @deprecated use symbolForLifeScienceID() instead
   * @param name of the string whose symbol you want to get
   * @throws NoSuchElementException if the string name is invalid.
   * @return The Symbol
   */
  static public Symbol symbolForName(String name)
  throws NoSuchElementException {
    String ls = "urn:lsid:biojava.org:symbol:"+name;
    LifeScienceIdentifier lsid = null;
    try {
      lsid = LifeScienceIdentifier.valueOf(ls);
    } catch (LifeScienceIdentifierParseException ex) {
      throw new BioError("Cannot construct LSID for "+name, ex);
    }
    Symbol s = (Symbol) lsidToSymbol.get(lsid);
    if(s == null) {
      throw new NoSuchElementException("Could not find symbol under the name " + lsid);
    }
    return s;
  }

  /**
   * Retreives the Symbol for the LSID
   * @param lsid the URN for the Symbol
   * @return a reference to the Symbol
   */
  static public Symbol symbolForLifeScienceID(LifeScienceIdentifier lsid){
    return (Symbol)lsidToSymbol.get(lsid);
  }

  /**
   * Register an alphabet by name.
   *
   * @param name  the name by which it can be retrieved
   * @param alphabet the Alphabet to store
   */
  static public void registerAlphabet(String name, Alphabet alphabet) {
    nameToAlphabet.put(name, alphabet);
    if(alphabet instanceof AbstractAlphabet){ //this might be needed for serialization
          ((AbstractAlphabet)alphabet).setRegistered(true);
    }
  }
  
  /**
   * Register and Alphabet by more than one name. This allows aliasing
   * of an alphabet with two or more names. It is equivalent to calling
   * <code>registerAlphabet(String name, Alphabet alphabet)</code> several
   * times.
   *
   * @since 1.4
   * @param names  the names by which it can be retrieved
   * @param alphabet the Alphabet to store
   */
  static public void registerAlphabet(String[] names, Alphabet alphabet){
      for(int i = 0; i < names.length; i++){
          registerAlphabet(names[i], alphabet);
      }
  }
  
  /**
   * A set of names under which Alphabets have been registered.
   * @return a <code>Set</code> of <code>Strings</code>
   */
  static public Set registrations(){
      return Collections.unmodifiableSet(nameToAlphabet.keySet());
  }

  /**
   * Has an Alphabet been registered by that name
   * @param name the name of the alphabet
   * @return true if it has or false otherwise
   */
  static public boolean registered(String name){
    return nameToAlphabet.containsKey(name);
  }

  /**
   * Get an iterator over all alphabets known.
   *
   * @return an Iterator over Alphabet objects
   */
  static public Iterator alphabets() {
    return Collections.unmodifiableCollection(nameToAlphabet.values()).iterator();
  }

  /**
   * <p>
   * Get the special `gap' Symbol.
   * </p>
   *
   * <p>
   * The gap symbol is a Symbol that has an empty alphabet of matches. As such
   *, ever alphabet contains gap, as there is no symbol that matches gap, so
   * there is no case where an alphabet doesn't contain a symbol that matches
   * gap.
   * </p>
   *
   * <p>
   * Gap can be thought of as an empty sub-space within the space of all
   * possible symbols. If you are working in a cross-product alphabet, you
   * should chose whether to use gap to represent 'no symbol', or a basis symbol
   * of the appropriate size built entirely of gaps to represent 'no symbol in
   * each of the slots'. Perhaps this could be explained better.
   * </p>
   *
   * @return the system-wide symbol that represents a gap
   */
  static public Symbol getGapSymbol() {
    return gapSymbol;
  }

  /**
   * <p>
   * Get the gap symbol appropriate to this list of alphabets.
   * </p>
   *
   * <p>
   * The gap symbol with have the same shape a the alphabet list. It will be as
   * long as the list, and if any of the alphabets in the list have a dimension
   * greater than 1, it will also insert the appropriate gap there.
   * </p>
   *
   * @param alphas  List of alphabets
   * @return the appropriate gap symbol for the alphabet list
   */
  static public Symbol getGapSymbol(List alphas) {
    SizeQueen sq = new SizeQueen(alphas);
    Symbol s = (Symbol) gapBySize.get(sq);

    if(s == null) {
      if(alphas.size() == 0) { // should never be needed
        s = gapSymbol;
      } else if(alphas.size() == 1) { // should never happen
        Alphabet a = (Alphabet) alphas.get(0);
        s = getGapSymbol(a.getAlphabets());
      } else {
        List symList = new ArrayList(alphas.size());
        for(Iterator i = alphas.iterator(); i.hasNext(); ) {
          Alphabet a = (Alphabet) i.next();
          symList.add(getGapSymbol(a.getAlphabets()));
        }
        try {
          s = new WellKnownGapSymbol(symList, sq);
        } catch (IllegalSymbolException ise) {
          throw new BioError(
            "Assertion Failure: Should be able to make gap basis", ise
          );
        }
      }
      gapBySize.put(sq, s);
    }

    return s;
  }
  
  

  /**
   * <p>
   * Generate a new AtomicSymbol instance with a name and Annotation.
   * </p>
   *
   * <p>
   * Use this method if you wish to create an AtomicSymbol instance. Initially it
   * will not be a member of any alphabet.
   * </p>
   *
   * @param name  the String returned by getName()
   * @param annotation the Annotation returned by getAnnotation()
   * @return a new AtomicSymbol instance
   */
  static public AtomicSymbol createSymbol(
    String name, Annotation annotation
  ) {
    AtomicSymbol as = new FundamentalAtomicSymbol(name, annotation);
    return as;
  }

  /**
   * <p>
   * Generate a new AtomicSymbol instance with a name and an Empty Annotation.
   * </p>
   *
   * <p>
   * Use this method if you wish to create an AtomicSymbol instance. Initially it
   * will not be a member of any alphabet.
   * </p>
   *
   * @param name  the String returned by getName()
   * @return a new AtomicSymbol instance
   */
  static public AtomicSymbol createSymbol(
      String name
      ) {
    AtomicSymbol as = new FundamentalAtomicSymbol(name, Annotation.EMPTY_ANNOTATION);
    return as;
  }

  /**
   * <p>
   * Generate a new AtomicSymbol instance with a token, name and Annotation.
   * </p>
   *
   * <p>
   * Use this method if you wish to create an AtomicSymbol instance. Initially it
   * will not be a member of any alphabet.
   * </p>
   *
   * @param token  the Char token returned by getToken() (ignpred as of BioJava 1.2)
   * @param name  the String returned by getName()
   * @param annotation the Annotation returned by getAnnotation()
   * @return a new AtomicSymbol instance
   * @deprecated Use the two-arg version of this method instead.
   */
  static public AtomicSymbol createSymbol(
    char token, String name, Annotation annotation
  ) {
    AtomicSymbol as = new FundamentalAtomicSymbol(name, annotation);
    return as;
  }

  /**
   * <p>
   * Generates a new Symbol instance that represents the tuple of Symbols in
   * symList.
   * </p>
   * 
   * <p>
   * This method is most useful for writing Alphabet implementations. It should
   * not be invoked by casual users. Use alphabet.getSymbol(List) instead.
   * </p>
   * @return a Symbol that encapsulates that List
   * @deprecated use the new version, without the token argument
   * @param annotation The annotation bundle for the symbol
   * @param token the Symbol's token [ignored since 1.2]
   * @param symList a list of Symbol objects
   * @param alpha the Alphabet that this Symbol will reside in
   * @throws org.biojava.bio.symbol.IllegalSymbolException If the Symbol cannot be made
   */
  static public Symbol createSymbol(
    char token, Annotation annotation,
    List symList, Alphabet alpha
  ) throws IllegalSymbolException {
      return createSymbol(annotation, symList, alpha);
  }

  static private Symbol readFromCache(List symList)
  {
    //System.out.println("Reading symbol: " + symList + " -> " + symListToSymbol.get(symList));
    return (Symbol) symListToSymbol.get(symList);
  }

  static private void writeToCache(List symList, Symbol sym)
  {
    //System.out.println("Writing symbol: " + symList + " -> " + sym);
    symListToSymbol.put(new ArrayList(symList), sym);
  }

  /**
   * <p>
   * Generates a new Symbol instance that represents the tuple of Symbols in
   * symList. This will attempt to return the same symbol for the same list.
   * </p>
   * 
   * <p>
   * This method is most useful for writing Alphabet implementations. It should
   * not be invoked by casual users. Use alphabet.getSymbol(List) instead.
   * </p>
   * @return a Symbol that encapsulates that List
   * @param annotation The annotation bundle for the Symbol
   * @param symList a list of Symbol objects
   * @param alpha the Alphabet that this Symbol will reside in
   * @throws org.biojava.bio.symbol.IllegalSymbolException If the Symbol cannot be made
   */
  static public Symbol createSymbol(
    Annotation annotation,
    List symList, Alphabet alpha)
          throws IllegalSymbolException
  {
    Symbol cs = readFromCache(symList);
    if(cs != null) {
      return cs;
    }

    Iterator i = symList.iterator();
    int basis = 0;
    int atomC = 0;
    int gaps = 0;
    while(i.hasNext()) {
      Symbol s = (Symbol) i.next();
      if(s instanceof BasisSymbol) {
        basis++;
        if(s instanceof AtomicSymbol) {
          atomC++;
        }
      } else {
        Alphabet matches = s.getMatches();
        if(matches instanceof FiniteAlphabet) {
          if(((FiniteAlphabet) matches).size() == 0) {
            gaps++;
          }
        }
      }
    }

    try {
      if(atomC == symList.size()) {
        Symbol sym = new SimpleAtomicSymbol(annotation, symList);
        writeToCache(symList, sym);
        return sym;
      } else if((gaps + basis) == symList.size()) {
        Symbol sym = new SimpleBasisSymbol(
                annotation,
                symList,
                new SimpleAlphabet(
                        expandMatches(alpha, symList, new ArrayList())));
        writeToCache(symList, sym);
        return sym;
      } else {
        Symbol sym = new SimpleSymbol(
                annotation,
                new SimpleAlphabet(
                        expandBasis(alpha, symList, new ArrayList())));
        writeToCache(symList,  sym);
        return sym;
      }
    } catch (IllegalSymbolException ise) {
      throw new IllegalSymbolException(
              ise,
              "Could not create a new symbol with: " +
              annotation + "\t" +
              symList + "\t" +
              alpha);
    }
  }

  /**
   * Expands a list of BasisSymbols into the set of AtomicSymbol instances
   * it matches.
   */
  private static Set expandBasis(Alphabet alpha, List symList, List built) {
    int indx = built.size();
    if(indx < symList.size()) {
      Symbol s = (Symbol) symList.get(indx);
      if(s instanceof AtomicSymbol) {
        built.add(s);
        return expandBasis(alpha, symList, built);
      } else {
        Set res = new HashSet();
        Iterator i = ((FiniteAlphabet) s.getMatches()).iterator();
        while(i.hasNext()) {
          AtomicSymbol as = (AtomicSymbol) i.next();
          List built2 = new ArrayList(built);
          built2.add(as);
          res.addAll(expandBasis(alpha, symList, built2));
        }
        return res;
      }
    } else {
      try {
        return Collections.singleton(alpha.getSymbol(built));
      } catch (IllegalSymbolException ise) {
        throw new BioError(
          "Assertion Failure: Should just have legal AtomicSymbol instances.", ise
        );
      }
    }
  }

  /**
   * <p>
   * Generates a new Symbol instance that represents the tuple of Symbols in
   * symList.
   * </p>
   * 
   * <p>
   * This method is most useful for writing Alphabet implementations. It should
   * not be invoked by users. Use alphabet.getSymbol(Set) instead.
   * </p>
   * @return a Symbol that encapsulates that List
   * @deprecated use the three-arg version of this method instead.
   * @param token the Symbol's token [ignored since 1.2]
   * @param annotation the Symbol's Annotation
   * @param symSet a Set of Symbol objects
   * @param alpha the Alphabet that this Symbol will reside in
   * @throws org.biojava.bio.symbol.IllegalSymbolException If the Symbol cannot be made
   */
  static public Symbol createSymbol(
    char token, Annotation annotation,
    Set symSet, Alphabet alpha
  ) throws IllegalSymbolException {
      return createSymbol(annotation, symSet, alpha);
  }

  /**
   * <p>
   * Generates a new Symbol instance that represents the tuple of Symbols in
   * symList.
   * </p>
   * 
   * <p>
   * This method is most useful for writing Alphabet implementations. It should
   * not be invoked by users. Use alphabet.getSymbol(Set) instead.
   * </p>
   * @return a Symbol that encapsulates that List
   * @param annotation the Symbol's Annotation
   * @param symSet a Set of Symbol objects
   * @param alpha the Alphabet that this Symbol will reside in
   * @throws org.biojava.bio.symbol.IllegalSymbolException If the Symbol cannot be made
   */
  static public Symbol createSymbol(
    Annotation annotation,
    Set symSet, Alphabet alpha
  ) throws IllegalSymbolException {
    if(symSet.size() == 0) {
      return getGapSymbol();
    }
    Set asSet = new HashSet();
    int len = -1;
    for(
      Iterator i = symSet.iterator();
      i.hasNext();
    ) {
      Symbol s = (Symbol) i.next();
      if(s instanceof AtomicSymbol) {
        AtomicSymbol as = (AtomicSymbol) s;
        int l = as.getSymbols().size();
        if(len == -1) {
          len = l;
        } else if(len != l) {
          throw new IllegalSymbolException(
            "Can't build ambiguity symbol as the symbols have inconsistent " +
            "length"
          );
        }
        asSet.add(as);
      } else {
        for(Iterator j = ((FiniteAlphabet) s.getMatches()).iterator();
          j.hasNext();
        ) {
          AtomicSymbol as = ( AtomicSymbol) j.next();
          int l = as.getSymbols().size();
          if(len == -1) {
            len = l;
          } else if(len != l) {
            throw new IllegalSymbolException(
              "Can't build ambiguity symbol as the symbols have inconsistent " +
              "length"
            );
          }
          asSet.add(as);
        }
      }
    }
    if(asSet.size() == 0) {
      return getGapSymbol();
    } else if(asSet.size() == 1) {
      return (Symbol) asSet.iterator().next();
    } else {
      if(len == 1) {
        return new SimpleBasisSymbol(
          annotation, new SimpleAlphabet(asSet)
        );
      } else {
        List fs = factorize(alpha, asSet);
        if(fs == null) {
          return new SimpleSymbol(
            annotation,
            new SimpleAlphabet(asSet)
          );
        } else {
          return new SimpleBasisSymbol(
            annotation,
            fs, new SimpleAlphabet(
              expandBasis(alpha, fs, new ArrayList())
            )
          );
        }
      }
    }
  }

  /**
   * Generates a new CrossProductAlphabet from the give name.
   *
   * @param name  the name to parse
   * @return the associated Alphabet
   */
  static public Alphabet generateCrossProductAlphaFromName(
    String name
  ) {
    if(!name.startsWith("(") || !name.endsWith(")")) {
      throw new BioError(
        "Can't parse " + name +
        " into a cross-product alphabet as it is not bracketed"
      );
    }

    name = name.substring(1, name.length()-1).trim();
    List aList = new ArrayList(); // the alphabets
    int i = 0;
    while(i < name.length()) {
      if(name.charAt(i) == '(') {
        int depth = 1;
        int j = i+1;
        while(j < name.length() && depth > 0) {
          char c = name.charAt(j);
          if(c == '(') {
            depth++;
          } else if(c == ')') {
            depth--;
          }
          j++;
        }
        if(depth == 0) {
          aList.add(alphabetForName(name.substring(i, j)));
          i = j;
        } else {
          throw new BioError(
            "Error parsing alphabet name: could not find matching bracket\n" +
            name.substring(i)
          );
        }
      } else {
        int j = name.indexOf(" x ", i);
        if(j < 0) {
          aList.add(alphabetForName(name.substring(i).trim()));
          i = name.length();
        } else {
          if(i != j){
            aList.add(alphabetForName(name.substring(i, j).trim()));
          }
          i = j + " x ".length();
        }
      }
    }

    return getCrossProductAlphabet(aList);
  }

  /**
   * <p>
   * Retrieve a CrossProductAlphabet instance over the alphabets in aList.
   * </p>
   * 
   * <p>
   * If all of the alphabets in aList implements FiniteAlphabet then the
   * method will return a FiniteAlphabet. Otherwise, it returns a non-finite
   * alphabet.
   * </p>
   * 
   * <p>
   * If you call this method twice with a list containing the same alphabets,
   * it will return the same alphabet. This promotes the re-use of alphabets
   * and helps to maintain the 'flyweight' principal for finite alphabet
   * symbols.
   * </p>
   * 
   * <p>
   * The resulting alphabet cpa will be retrievable via
   * AlphabetManager.alphabetForName(cpa.getName())
   * </p>
   * @param aList a list of Alphabet objects
   * @return a CrossProductAlphabet that is over the alphabets in aList
   */
  static public Alphabet getCrossProductAlphabet(List aList) {
    return getCrossProductAlphabet(aList, (Alphabet) null);
  }

  
  /**
   * Attempts to create a cross product alphabet and register it under a name.
   * @param aList A list of alphabets
   * @param name The name which the new alphabet will be registered under.
   * @throws org.biojava.bio.symbol.IllegalAlphabetException If the Alphabet cannot be made or a different 
   * alphabet is already registed under this name.
   * @return The CrossProductAlphabet
   */
  static public Alphabet getCrossProductAlphabet(List aList, String name)
  throws IllegalAlphabetException {
    Alphabet currentAlpha = (Alphabet) nameToAlphabet.get(name);
    if(currentAlpha != null) {
      if(currentAlpha.getAlphabets().equals(aList)) {
        return currentAlpha;
      } else {
        throw new IllegalAlphabetException(name + " already registered");
      }
    } else {
      Alphabet alpha = getCrossProductAlphabet(aList);
      registerAlphabet(name, alpha);
      return alpha;
    }
  }

  /**
   * <p>
   * Retrieve a CrossProductAlphabet instance over the alphabets in aList.
   * </p>
   *
   * <p>
   * This method is most usefull for implementors of cross-product alphabets,
   * allowing them to safely build the matches alphabets for ambiguity symbols.
   * </p>
   *
   * <p>
   * If all of the alphabets in aList implements FiniteAlphabet then the
   * method will return a FiniteAlphabet. Otherwise, it returns a non-finite
   * alphabet.
   * </p>
   *
   * <p>
   * If you call this method twice with a list containing the same alphabets,
   * it will return the same alphabet. This promotes the re-use of alphabets
   * and helps to maintain the 'flyweight' principal for finite alphabet
   * symbols.
   * </p>
   *
   * <p>
   * The resulting alphabet cpa will be retrievable via
   * AlphabetManager.alphabetForName(cpa.getName())
   * </p>
   *
   * @param aList a list of Alphabet objects
   * @param parent a parent alphabet
   * @return a CrossProductAlphabet that is over the alphabets in aList
   */
  static public Alphabet getCrossProductAlphabet(
    List aList, Alphabet parent
  ) {
    if(aList.size() == 0) {
      return Alphabet.EMPTY_ALPHABET;
    }

    // This trap means that the `product' operator can be
    // safely applied to a single alphabet.

    if (aList.size() == 1)
        return (Alphabet) aList.get(0);

    if(crossProductAlphabets == null) {
      crossProductAlphabets = new HashMap();
    }

    Alphabet cpa = (Alphabet) crossProductAlphabets.get(aList);

    int size = 1;
    if(cpa == null) {
      for(Iterator i = aList.iterator(); i.hasNext(); ) {
        Alphabet aa = (Alphabet) i.next();
        if(! (aa instanceof FiniteAlphabet) ) {
          cpa =  new InfiniteCrossProductAlphabet(aList);
          break;
        }
        if(size <= 1000) {
          size *= ((FiniteAlphabet) aa).size();
        }
      }
      if(cpa == null) {
        try {
          if(size > 0 && size < 1000) {
            cpa = new SimpleCrossProductAlphabet(aList, parent);
          } else {
            cpa = new SparseCrossProductAlphabet(aList);
          }
        } catch (IllegalAlphabetException iae) {
          throw new BioError(
            "Could not create SimpleCrossProductAlphabet for " + aList +
            " even though we should be able to. No idea what is wrong."
          );
        }
      }
      crossProductAlphabets.put(new ArrayList(aList), cpa);
      registerAlphabet(cpa.getName(), cpa);
    }

    return cpa;
  }

  private static Set expandMatches(Alphabet parent, List symList, List built) {
    int indx = built.size();
    if(indx < symList.size()) {
      Symbol bs = (Symbol) symList.get(indx);
      if(bs instanceof AtomicSymbol) {
        built.add(bs);
        return expandMatches(parent, symList, built);
      } else {
        Set syms = new HashSet();
        Iterator i = ((FiniteAlphabet) bs.getMatches()).iterator();
        while(i.hasNext()) {
          List built2 = new ArrayList(built);
          built2.add((AtomicSymbol) i.next());
          syms.addAll(expandMatches(parent, symList, built2));
        }
        return syms;
      }
    } else {
      try {
        Symbol s = parent.getSymbol(built);
        if(s instanceof AtomicSymbol) {
          return Collections.singleton((AtomicSymbol) s);
        } else {
          Set syms = new HashSet();
          for(Iterator i = ((FiniteAlphabet) s.getMatches()).iterator(); i.hasNext(); ) {
            syms.add((AtomicSymbol) i.next());
          }
          return syms;
        }
      } catch (IllegalSymbolException ise) {
        throw new BioError("Assertion Failure: Couldn't create symbol.", ise);
      }
    }
  }

  /**
   * <p>
   * Return a list of BasisSymbol instances that uniquely sum up all
   * AtomicSymbol
   * instances in symSet. If the symbol can't be represented by a single list of
   * BasisSymbol instances, return null.
   * </p>
   * 
   * <p>
   * This method is most useful for implementers of Alphabet and Symbol. It
   * probably should not be invoked by users.
   * </p>
   * @return a List of BasisSymbols
   * @param symSet the Set of AtomicSymbol instances
   * @param alpha the Alphabet instance that the Symbols are from
   * @throws org.biojava.bio.symbol.IllegalSymbolException In practice it should not. If it does it probably
   * indicates a subtle bug somewhere in AlphabetManager
   */
  public static List factorize(Alphabet alpha, Set symSet)
  throws IllegalSymbolException {
    List alphas = alpha.getAlphabets();
    List facts = new ArrayList();
    int size = symSet.size();
    Set syms = new HashSet();
    for(int col = 0; col < alphas.size(); col++) {
      Alphabet a = (Alphabet) alphas.get(col);
      for(Iterator i = symSet.iterator(); i.hasNext(); ) {
        syms.add(
          (AtomicSymbol) ((AtomicSymbol)
          i.next()).getSymbols().get(col)
        );
      }
      int s = syms.size();
      if( (size % s) != 0 ) {
        return null;
      }
      size /= s;
      facts.add(a.getAmbiguity(syms));
      syms.clear();
    }
    if(size != 1) {
      return null;
    }
    return facts;
  }




    /**
     * Load additional Alphabets, defined in XML format, into the AlphabetManager's registry.
     * These can the be retrieved by calling <code>alphabetForName</code>.
     *
     * @param is an <code>InputSource</code> encapsulating the document to be parsed
     * @throws IOException if there is an error accessing the stream
     * @throws SAXException if there is an error while parsing the document
     * @throws BioException if a problem occurs when creating the new Alphabets.
     * @since 1.3
     */

    public static void loadAlphabets(InputSource is)
        throws SAXException, IOException, BioException
    {
        try {
            SAXParserFactory spf = SAXParserFactory.newInstance();
            spf.setNamespaceAware(true);
            XMLReader parser = spf.newSAXParser().getXMLReader();
            parser.setContentHandler(new SAX2StAXAdaptor(new AlphabetManagerHandler()));
            parser.parse(is);
        } catch (ParserConfigurationException ex) {
            throw new BioException( "Unable to create XML parser", ex);
        }
    }

    /**
     * StAX handler for the alphabetManager element
     */

    private static class AlphabetManagerHandler extends StAXContentHandlerBase {
        public void startElement(String nsURI,
                                             String localName,
                                             String qName,
                                             Attributes attrs,
                                             DelegationManager dm)
             throws SAXException
         {
             if (localName.equals("alphabetManager")) {
                 // ignore
             } else if (localName.equals("symbol")) {
                 String name = attrs.getValue("name");
                 dm.delegate(new SymbolHandler(name));
             } else if (localName.equals("alphabet")) {
                 String name = attrs.getValue("name");
                 String parent = attrs.getValue("parent");
                 FiniteAlphabet parentAlpha = null;
                 if (parent != null && parent.length() > 0) {
                     parentAlpha = (FiniteAlphabet) nameToAlphabet.get(parent);
                 }
                 dm.delegate(new AlphabetHandler(name, parentAlpha));
             } else {
                 throw new SAXException(
                         "Unknown element in alphabetManager: " +
                         localName);
             }
         }

         public void endElement(String nsURI,
                                String localName,
                                String qName,
                                StAXContentHandler delegate)
            throws SAXException
         {
             if (delegate instanceof SymbolHandler) {
                 SymbolHandler sh = (SymbolHandler) delegate;
                 //String name = sh.getName();
                 LifeScienceIdentifier lsid = sh.getLSID();
                 Symbol symbol = sh.getSymbol();
                 if (lsidToSymbol.containsKey(lsid)) {
                     throw new SAXException(
                     "There is already a top-level symbol named "
                     + lsid);
                 }
                 lsidToSymbol.put(lsid, symbol);
             } else if (delegate instanceof AlphabetHandler) {
                 AlphabetHandler ah = (AlphabetHandler) delegate;
                 String name = ah.getName();
                 FiniteAlphabet alpha = ah.getAlphabet();
                 registerAlphabet(name, alpha);
             }
         }

         private class SymbolHandler extends StAXContentHandlerBase {
             private String name;
             private LifeScienceIdentifier lsid;
             private Symbol symbol;
             private Annotation annotation = new SmallAnnotation();

             public SymbolHandler(String id) {
                try {
                  lsid = LifeScienceIdentifier.valueOf(id);
                  name = lsid.getObjectId();
                } catch (LifeScienceIdentifierParseException ex) {
                  throw new BioError("Malformed LSID - "+name, ex);
                }
             }

             public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                  throws SAXException
             {
                 if (localName.equals("symbol")) {
                     // ignore
                 } else if (localName.equals("description")) {
                     dm.delegate(new StringElementHandlerBase() {
                         protected void setStringValue(String s) {
                             try {
                                 annotation.setProperty("description", s);
                             } catch (ChangeVetoException ex) {
                                 throw new BioError( "Assertion failure: veto while modifying new Annotation", ex);
                             }
                         }
                     } );
                 } else {
                     throw new SAXException("Unknown element in symbol: " + localName);
                 }
             }

             public void endTree() {
                 symbol = new WellKnownAtomicSymbol(
                    new FundamentalAtomicSymbol(
                        name,
                        annotation
                    ),
                    lsid
                  );
             }

             Symbol getSymbol() {
                 return symbol;
             }

             String getName() {
                 return name;
             }

             LifeScienceIdentifier getLSID(){
               return lsid;
             }
         }

         private class AlphabetHandler extends StAXContentHandlerBase {
             private String name;
             //private Map localSymbols;
             private WellKnownAlphabet alpha;
             private ImmutableWellKnownAlphabetWrapper alphaWrapper;

             String getName() {
                 return name;
             }

             FiniteAlphabet getAlphabet() {
                 return alphaWrapper;
             }

             public void endTree() {
                 alpha.addChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN);
             }

             public AlphabetHandler(String name, FiniteAlphabet parent) {
                 this.name = name;
                 //localSymbols = new OverlayMap(nameToSymbol);
                 alpha = new WellKnownAlphabet();
                 alpha.setName(name);
                 alphaWrapper = new ImmutableWellKnownAlphabetWrapper(alpha);
                 if (parent != null) {
                     for (Iterator i = parent.iterator(); i.hasNext(); ) {
                         WellKnownAtomicSymbol sym =
                             (WellKnownAtomicSymbol) i.next();
                         try {
                             alpha.addSymbol(sym);
                         } catch (Exception ex) {
                             throw new BioError(
                               "Couldn't initialize alphabet from parent", ex);
                         }
                         lsidToSymbol.put(sym.getIdentifier(), sym);
                     }
                 }
             }

             public void startElement(String nsURI,
                                     String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                  throws SAXException
             {
                 if (localName.equals("alphabet")) {
                     // ignore
                 } else if (localName.equals("symbol")) {
                     String name = attrs.getValue("name");
                     dm.delegate(new SymbolHandler(name));
                 } else if (localName.equals("symbolref")) {
                     String name = attrs.getValue("name");
                    LifeScienceIdentifier lsid = null;
                    try {
                      lsid =
                          LifeScienceIdentifier.valueOf(name);
                    } catch (LifeScienceIdentifierParseException ex) {
                      throw new SAXException("Couldn't form a LSID from "+name);
                    }
                     Symbol sym = (Symbol) lsidToSymbol.get(lsid);
                     if (sym == null) {
                         throw new SAXException(
                           "Reference to non-existent symbol " + name);
                     }
                     addSymbol(sym);
                 } else if (localName.equals("characterTokenization")) {
                     String name = attrs.getValue("name");
                     boolean caseSensitive = "true".equals(attrs.getValue("caseSensitive"));
                     dm.delegate(new CharacterTokenizationHandler(name, alphaWrapper, lsidToSymbol, caseSensitive));
                 } else if (localName.equals("description")) {
                     dm.delegate(new StringElementHandlerBase() {
                         protected void setStringValue(String s) {
                             try {
                                 alpha.getAnnotation().setProperty("description", s);
                             } catch (ChangeVetoException ex) {
                                 throw new BioError( "Assertion failure: veto while modifying new Annotation", ex);
                             }
                         }
                     } );
                 } else {
                     throw new SAXException("Unknown element in alphabetl: " + localName);
                 }
             }

             public void endElement(String nsURI,
                                                String localName,
                                    String qName,
                                    StAXContentHandler delegate)
                  throws SAXException
             {
                 if (delegate instanceof SymbolHandler) {
                     SymbolHandler sh = (SymbolHandler) delegate;
                     //String name = sh.getName();
                     Symbol symbol = sh.getSymbol();
                     LifeScienceIdentifier lsid = sh.getLSID();
                     lsidToSymbol.put(lsid, symbol);
                     addSymbol(symbol);
                 } else if (delegate instanceof CharacterTokenizationHandler) {
                     CharacterTokenizationHandler cth = (CharacterTokenizationHandler) delegate;
                     String name = cth.getName();
                     SymbolTokenization toke = cth.getTokenization();
                     alpha.putTokenization(name, toke);
                 }
             }

             private void addSymbol(Symbol sym)
                 throws SAXException
             {
                 try {
                     alpha.addSymbol(sym);
                 } catch (ChangeVetoException cve) {
                     throw new BioError( "Assertion failure: veto while modifying new Alphabet", cve);
                 } catch (IllegalSymbolException ex) {
                     throw new SAXException("IllegalSymbolException adding symbol to alphabet");
                 }
             }
         }

         private class CharacterTokenizationHandler extends StAXContentHandlerBase {
             private String name;
             private Map localSymbols;
             private SymbolTokenization toke;
             private boolean isAlternate;

             String getName() {
                 return name;
             }

             SymbolTokenization getTokenization() {
                 return toke;
             }

             public CharacterTokenizationHandler(String name,
                                                 FiniteAlphabet alpha,
                                                 Map localSymbols,
                                                 boolean caseSensitive)
             {

                 this.name = name;
                 this.localSymbols = new HashMap();
                 for (Iterator i = alpha.iterator(); i.hasNext(); ) {
                     WellKnownAtomicSymbol sym = (WellKnownAtomicSymbol) i.next();
                     this.localSymbols.put(sym.getIdentifier(), sym);
                 }
                 if(name.indexOf("alternate")==0) {
                     toke = new AlternateTokenization(alpha, caseSensitive);
                     isAlternate = true;
                 } else
                     toke = new CharacterTokenization(alpha, caseSensitive);
             }

             public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                  throws SAXException
             {
                 if (localName.equals("characterTokenization")) {
                     // ignore
                 } else if (localName.equals("atomicMapping")) {
                     dm.delegate(new MappingHandler(true));
                 } else if (localName.equals("ambiguityMapping")) {
                     dm.delegate(new MappingHandler(false));
                 } else if (localName.equals("gapSymbolMapping")) {
                     dm.delegate(new MappingHandler(false, true));
                 } else {
                     throw new SAXException("Unknown element in characterTokenization: " + localName);
                 }
             }

             private class MappingHandler extends StAXContentHandlerBase {
                 public MappingHandler(boolean isAtomic, boolean isPureGap) {
                   this.isAtomic = isAtomic;
                   this.isPureGap = isPureGap;
                 }

                 public MappingHandler(boolean isAtomic) {
                     this(isAtomic, false);
                 }

                 boolean isAtomic;
                 boolean isPureGap;
                 Set symbols = new HashSet();
                 char c = '\0';
                 String str = "";
                 int level = 0;

                 public void startElement(String nsURI,
                                          String localName,
                                          String qName,
                                          Attributes attrs,
                                          DelegationManager dm)
                     throws SAXException
                 {
                     if (level == 0) {
                         c = attrs.getValue("token").charAt(0);
                         if(isAlternate)
                             str = attrs.getValue("token");
                     } else {
                         if (localName.equals("symbolref")) {
                             String name = attrs.getValue("name");
                             LifeScienceIdentifier lsid = null;
                             try {
                               lsid = LifeScienceIdentifier.valueOf(name);
                             } catch (LifeScienceIdentifierParseException ex) {
                               throw new SAXException("Cannot for LSID from " + name);
                             }
                             Symbol sym = (Symbol) localSymbols.get(lsid);
                             if (sym == null) {
                                 throw new SAXException("Reference to non-existent symbol " + name);
                             }
                             symbols.add(sym);
                         } else {
                             throw new SAXException("Unknown element in mapping: " + localName);
                         }
                     }
                     ++level;
                 }

                 public void endElement(String nsURI,
                                        String localName,
                                        String qName,
                                        StAXContentHandler delegate)
                     throws SAXException
                 {
                     --level;
                 }

                 public void endTree()
                     throws SAXException
                 {
                     Symbol ambiSym;
                     if(isPureGap) {
                         ambiSym = getGapSymbol();
                     } else {
                         try {
                             ambiSym = toke.getAlphabet().getAmbiguity(symbols);
                         } catch (IllegalSymbolException ex) {
                             throw (SAXException)
                                     new SAXException("IllegalSymbolException binding mapping for " + c).initCause(ex);
                         }
                     }
                     if(isAlternate)
                        ((AlternateTokenization)toke).bindSymbol(ambiSym, str);
                     else
                        ((CharacterTokenization)toke).bindSymbol(ambiSym, c);
                 }
             }
         }
    }

    private static class WellKnownTokenizationWrapper
        extends Unchangeable
        implements SymbolTokenization, Serializable
    {
        private String name;
        private Alphabet alphabet;
        private SymbolTokenization toke;

        WellKnownTokenizationWrapper(Alphabet alpha, SymbolTokenization toke, String name) {
            super();
            this.alphabet = alpha;
            this.name = name;
            this.toke = toke;
        }

        public Alphabet getAlphabet() {
            return alphabet;
        }

        public TokenType getTokenType() {
            return toke.getTokenType();
        }

        public StreamParser parseStream(SeqIOListener listener) {
            return toke.parseStream(listener);
        }

        public Symbol parseToken(String s)
            throws IllegalSymbolException
        {
            return toke.parseToken(s);
        }

        public String tokenizeSymbol(Symbol s)
            throws IllegalSymbolException
        {
            return toke.tokenizeSymbol(s);
        }

        public String tokenizeSymbolList(SymbolList sl)
            throws IllegalAlphabetException, IllegalSymbolException
        {
            return toke.tokenizeSymbolList(sl);
        }

        public Annotation getAnnotation() {
            return toke.getAnnotation();
        }

        public Object writeReplace() {
            return new OPH(getAlphabet().getName(), name);
        }

        private static class OPH implements Serializable {
            private String alphaName;
            private String name;

            OPH(String alphaName, String name) {
                this.alphaName = alphaName;
                this.name = name;
            }

            private Object readResolve() throws ObjectStreamException {
                try {
                    Alphabet alphabet = alphabetForName(alphaName);
                    return alphabet.getTokenization(name);
                } catch (Exception ex) {
                    throw new InvalidObjectException("Couldn't resolve tokenization " + name + " in alphabet " + alphaName);
                }
            }
        }
    }

    /**
     * An alphabet contained WellKnownSymbols
     */

    private static class WellKnownAlphabet
        extends SimpleAlphabet
    {
        public WellKnownAlphabet() {
            super();
        }

        public WellKnownAlphabet(Set s) {
            super(s);
        }

        protected Symbol getAmbiguityImpl(Set s)
            throws IllegalSymbolException
        {
            return getWellKnownAmbiguitySymbol(s);
        }
    }

    /**
     * A wrapper which makes an Alphabet unchangable, and also fixes serialization
     */

    private static class ImmutableWellKnownAlphabetWrapper
        extends Unchangeable
        implements FiniteAlphabet, Serializable
    {
        private FiniteAlphabet alpha;
        private Map tokenizationsByName = new HashMap();

        public ImmutableWellKnownAlphabetWrapper(FiniteAlphabet alpha) {
            super();
            this.alpha = alpha;
        }

        private Object writeReplace() {
            return new OPH(getName());
        }

        public SymbolTokenization getTokenization(String name)
            throws BioException
        {
            SymbolTokenization toke = (SymbolTokenization) tokenizationsByName.get(name);
            if (toke == null) {
                if ("name".equals(name)) {
                    toke = new NameTokenization(this);
                } else {
                    toke = new WellKnownTokenizationWrapper(this, alpha.getTokenization(name), name);
                }
                tokenizationsByName.put(name, toke);
            }
            return toke;
        }

        /**
         * Placeholder for a WellKnownAlphabet in a serialized
         * object stream.
         */

         private static class OPH implements Serializable {
             private String name;

             public OPH(String name) {
                 this.name = name;
             }

             private Object readResolve() throws ObjectStreamException {
                 try {
                     Alphabet a = AlphabetManager.alphabetForName(name);
                     return a;
                 } catch (NoSuchElementException ex) {
                     throw new InvalidObjectException("Couldn't resolve alphabet " + name);
                 }
             }
         }

        public boolean contains(Symbol s) {
            return alpha.contains(s);
        }

        public List getAlphabets() {
            return Collections.singletonList(this);
        }

        public Symbol getAmbiguity(Set s)
            throws IllegalSymbolException
        {
            return alpha.getAmbiguity(s);
        }

        public Symbol getGapSymbol() {
            return alpha.getGapSymbol();
        }

        public String getName() {
            return alpha.getName();
        }

        public Symbol getSymbol(List l)
            throws IllegalSymbolException
        {
            return alpha.getSymbol(l);
        }

        public void validate(Symbol s)
            throws IllegalSymbolException
        {
                alpha.validate(s);
        }

        public void addSymbol(Symbol s)
            throws ChangeVetoException
        {
            throw new ChangeVetoException("Can't add symbols to Well Known Alphabets");
        }

        public void removeSymbol(Symbol s)
            throws ChangeVetoException
        {
            throw new ChangeVetoException("Can't remove symbols from Well Known Alphabets");
        }

        public Iterator iterator() {
            return  alpha.iterator();
        }

        public int size() {
            return alpha.size();
        }

        public Annotation getAnnotation() {
            return alpha.getAnnotation();
        }
    }

    
    /**
     * A well-known gap. Resolved in serialized data
     */
    private static class WellKnownGapSymbol extends AbstractSimpleBasisSymbol implements Serializable{
        private SizeQueen sq;
        public WellKnownGapSymbol(List symList, SizeQueen sq) throws IllegalSymbolException{
            super(Annotation.EMPTY_ANNOTATION,
            symList,
            Alphabet.EMPTY_ALPHABET);
            this.sq = sq;
        }
        
        private Object readResolve() throws ObjectStreamException{
            //System.out.println("ping!!");
            return AlphabetManager.getGapSymbol(sq.getAlphas());
        }
    }
    /**
     * A well-known symbol.  Replaced by a placeholder in
     * serialized data.
     */

    private static class WellKnownAtomicSymbol
        extends WellKnownBasisSymbol
        implements AtomicSymbol, Identifiable {

        LifeScienceIdentifier lsid;

        WellKnownAtomicSymbol(AtomicSymbol symbol, LifeScienceIdentifier lsid) {
            super(symbol);
            this.lsid = lsid;
        }

        public LifeScienceIdentifier getIdentifier(){
          return lsid;
        }

        public Alphabet getMatches() {
            return new SingletonAlphabet(this);
        }

        private Object writeReplace() {
            return new WellKnownAtomicSymbol.OPH(getIdentifier());
        }

        /**
         * Object Place Holder
         */
        private static class OPH implements Serializable {
            private LifeScienceIdentifier name;

            public OPH(LifeScienceIdentifier name) {
                this.name = name;
            }

            private Object readResolve() throws ObjectStreamException {
                try {
                    return symbolForLifeScienceID(name);
                } catch (NoSuchElementException ex) {
                    throw new InvalidObjectException(
                        "Couldn't resolve symbol:" + name
                    );
                }
            }
        }
    }

    private static class WellKnownBasisSymbol
            extends Unchangeable
            implements BasisSymbol, Serializable
    {
        protected BasisSymbol symbol;
        private Set matches;

        WellKnownBasisSymbol(BasisSymbol symbol) {
            super();
            symbol.addChangeListener(ChangeListener.ALWAYS_VETO, ChangeType.UNKNOWN); // Immutable
            this.symbol = symbol;
            this.matches = new HashSet();
            for (Iterator i = ((FiniteAlphabet) symbol.getMatches()).iterator(); i.hasNext(); ) {
                matches.add(i.next());
            }
        }

        Symbol getSymbol() {
            return symbol;
        }

        public int hashCode() {
            return symbol.hashCode();
        }

        public boolean equals(Object o) {
            if (o instanceof WellKnownBasisSymbol) {
                return symbol.equals(((WellKnownBasisSymbol) o).getSymbol());
            } else {
                return false;
            }
        }

        public String getName() {
            return symbol.getName();
        }

        public Alphabet getMatches() {
            return symbol.getMatches();
        }

        public List getSymbols() {
            return Collections.singletonList(this);
        }

        public Annotation getAnnotation() {
            return symbol.getAnnotation();
        }

        private Object writeReplace() {
            return new OPH(matches);
        }


        private static class OPH implements Serializable {
            private Set matches;

            public OPH(Set matches) {
                OPH.this.matches = matches;
            }

            private Object readResolve() /* throws ObjectStreamException */ {
                return getWellKnownAmbiguitySymbol(matches);
            }
        }
    }

  /**
   * <p>
   * The class representing the Gap symbol.
   * </p>
   *
   * <p>
   * The gap is quite special. It is an ambiguity symbol with an empty alphabet.
   * This means that it notionaly represents an unfilled slot in a sequence.
   * It should be a singleton, hence the
   * placement in AlphabetManager and also the method normalize.
   * </p>
   *
   * @author Matthew Pocock
   */
  private static class GapSymbol
    extends
      Unchangeable
    implements
      Symbol,
      Serializable
  {
      public GapSymbol() {
      }

      public String getName() {
          return "gap";
      }

      public char getToken() {
          return '-';
      }

      public Annotation getAnnotation() {
          return Annotation.EMPTY_ANNOTATION;
      }

      public Alphabet getMatches() {
          return Alphabet.EMPTY_ALPHABET;
      }
      
      
       private Object readResolve() throws ObjectStreamException {
           return AlphabetManager.getGapSymbol();
       }
  }
  

  /**
   * Get an indexer for a specified alphabet.
   *
   * @param alpha The alphabet to index
   * @return an AlphabetIndex instance
   * @since 1.1
   */

  /**
   * Get an indexer for a specified alphabet.
   *
   * @param alpha The alphabet to index
   * @return an AlphabetIndex instance
   * @since 1.1
   */
  public static AlphabetIndex getAlphabetIndex(
    FiniteAlphabet alpha
  ) {
    final int generateIndexSize = 160;
    AlphabetIndex ai = (AlphabetIndex) alphabetToIndex.get(alpha);
    if(ai == null) {
      int size = alpha.size();
      if(size <= generateIndexSize) {
        ai = new LinearAlphabetIndex(alpha);
      } else {
        if(alpha.getAlphabets().size() > 1) {
          ai = new CrossProductAlphabetIndex(alpha);
        } else {
          ai = new HashedAlphabetIndex(alpha);
        }
      }
      alphabetToIndex.put(alpha, ai);
    }
    return ai;
  }

  /**
   * Get an indexer for an array of symbols.
   *
   * @param syms the Symbols to index in that order
   * @return an AlphabetIndex instance
   * @since 1.1
   */
  public static AlphabetIndex getAlphabetIndex (
    Symbol[] syms
  ) throws IllegalSymbolException, BioException {
    return new LinearAlphabetIndex(syms);
  }

  private static final class SizeQueen extends AbstractList implements Serializable{
    private final List alphas;

    public SizeQueen(List alphas) {
      this.alphas = alphas;
    }

    public int size() {
      return alphas.size();
    }
    
    public List getAlphas(){
        return this.alphas;
    }

    public Object get(int pos) {
      Alphabet a = (Alphabet) alphas.get(pos);
      List al = a.getAlphabets();
      int size = al.size();
      if(size > 1) {
        return new SizeQueen(al);
      } else {
        return new Integer(size);
      }
    }
  }

  private static Symbol getWellKnownAmbiguitySymbol(Set s) {
      Symbol sym = (Symbol) ambiguitySymbols.get(s);
      if (sym == null) {
          SimpleAlphabet matchAlpha = new WellKnownAlphabet(s);
          sym = new WellKnownBasisSymbol(new SimpleBasisSymbol(Annotation.EMPTY_ANNOTATION, matchAlpha));
          ambiguitySymbols.put(new HashSet(s), sym);
      }
      return sym;
  }
}
