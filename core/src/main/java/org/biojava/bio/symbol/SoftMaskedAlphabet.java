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

import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.StreamParser;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IntegerAlphabet.IntegerSymbol;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ListTools;
import org.biojava.utils.Unchangeable;

/**
 * Soft masking is usually displayed by making the masked regions somehow
 * different from the non masked regions. Typically the masked regions are
 * lower case but other schemes could be invented. For example a softmasked
 * DNA sequence may look like this:<pre>
 *
 * >DNA_sequence
 * ATGGACGCTAGCATggtggtggtggtggtggtggtGCATAGCGAGCAAGTGGAGCGT
 *
 * </pre>
 * Where the lowercase regions are masked by low complexity.
 * <p>
 * <code>SoftMaskedAlphabet</code>s come with <code>SymbolTokenizers</code>
 * that understand how to read and write the softmasking. The interpretation
 * of what constitutes a masked region is governed by an implementation of
 * a <code>MaskingDetector</code>. The <code>DEFAULT</code> field of the
 * <code>MaskingDetector</code> interface defines lower case tokens as masked.

 * <p> Copyright (c) 2004 Novartis Institute for Tropical Diseases</p>
 * @author Mark Schreiber
 * @version 1.0
 */

public final class SoftMaskedAlphabet
    extends Unchangeable implements FiniteAlphabet{

  //used to indicate masking. 0 indicates no mask 1 indicates mask.
  private IntegerAlphabet.SubIntegerAlphabet binary;
  private FiniteAlphabet alpha;
  private String name;
  private FiniteAlphabet delegateAlpha;
  private MaskingDetector maskingDetector;

  private SoftMaskedAlphabet(FiniteAlphabet alpha, String name)
      throws IllegalAlphabetException{
    this.alpha = alpha;
    binary = IntegerAlphabet.getSubAlphabet(0,1);
    this.name = name;
    delegateAlpha = (FiniteAlphabet)AlphabetManager.getCrossProductAlphabet(
        new ListTools.Doublet(alpha, binary));
  }

  /**
   * Generates a soft masked Alphabet where lowercase tokens are assumed to be
   * soft masked.
   * @param alphaToMask for example the DNA alphabet.
   * @throws IllegalAlphabetException if it cannot be constructed
   * @return a reference to a singleton <code>SoftMaskedAlphabet</code>.
   */
  public static SoftMaskedAlphabet getInstance(FiniteAlphabet alphaToMask)
      throws IllegalAlphabetException {
    return getInstance(alphaToMask, MaskingDetector.DEFAULT);
  }

  /**
   * Creates a compound alphabet that is a hybrid of the alphabet that is to
   * be soft masked and a binary alphabet that indicates if any
   * <code>Symbol</code> is soft masked or not.
   *
   * @param alphaToMask for example the DNA alphabet.
   * @param maskingDetector to define masking behaivour
   * @throws IllegalAlphabetException if it cannot be constructed
   * @return a reference to a singleton <code>SoftMaskedAlphabet</code>.
   */
  public static SoftMaskedAlphabet getInstance(FiniteAlphabet alphaToMask,
                                               MaskingDetector maskingDetector)
      throws IllegalAlphabetException{
    String lookup = "Softmasked {"+alphaToMask.getName()+"}";
    if(AlphabetManager.registered(lookup)){
      return (SoftMaskedAlphabet)AlphabetManager.alphabetForName(lookup);
    }

    SoftMaskedAlphabet sma = new SoftMaskedAlphabet(alphaToMask, lookup);
    AlphabetManager.registerAlphabet(sma.getName(), sma);

    sma.maskingDetector = maskingDetector;
    return sma;
  }

  /**
   * Gets the <CODE>Alphabet</CODE> upon which masking is being applied
   * @return A <CODE>FiniteAlphabet</CODE>
   */
  public FiniteAlphabet getMaskedAlphabet(){
    return alpha;
  }

  /**
   * The compound alpha that holds the symbols used by this wrapper
   * @return a <code>FiniteAlphabet</code>
   */
  protected FiniteAlphabet getDelegate(){
    return delegateAlpha;
  }

  /**
   * The SoftMaskedAlphabet has no annotation
   * @return Annotation.EMPTY_ANNOTATION
   */
  public Annotation getAnnotation(){
    return Annotation.EMPTY_ANNOTATION;
  }

  /**
   * The name of the Alphabet
   * @return a <code>String</code> in the form of
   * <code>"Softmasked {"+alphaToMask.getName()+"}"</code>
   */
  public String getName(){
    return name;
  }

  /**
   * Gets the components of the <code>Alphabet</code>.
   * @return a <code>List</code> with two members, the first is the wrapped
   * <code>Alphabet</code> the second is the binary
   * <code>SubIntegerAlphabet</code>.
   */
  public List getAlphabets(){
    return new ListTools.Doublet(alpha, binary);
  }

  
  /**
   * Gets the compound symbol composed of the <code>Symbols</code> in the List.
   * The <code>Symbols</code> in the <code>List</code> must be from <code>alpha</code>
   * (defined in the constructor) and <code>SUBINTEGER[0..1]</code>
   * @return A <code>Symbol</code> from this alphabet.
   * @throws IllegalSymbolException if <code>l</code> is not as expected (see above)
   * @param l a <code>List</code> of <code>Symbols</code>
   */
  public Symbol getSymbol(List l) throws IllegalSymbolException {
    return delegateAlpha.getSymbol(l);
  }
  
  /**
   * This is not supported. Ambiguity should be handled at the level of the 
   * wrapped Alphabet. Use <code>getSymbol(List l)</code> instead and provide
   * it with an ambigutiy and a masking symbol.
   * @param s a <code>Set</code> of <code>Symbols</code>
   * @see #getSymbol(List l)
   * @throws UnsupportedOperationException
   */
  public Symbol getAmbiguity(Set s) throws UnsupportedOperationException {
    throw new UnsupportedOperationException(
        "Ambiguity should be handled at the level of the wrapped Alphabet");
  }

  public Symbol getGapSymbol(){
    return AlphabetManager.getGapSymbol(new ListTools.Doublet(alpha, binary));
  }

  public boolean contains(Symbol s){
    return delegateAlpha.contains(s);
  }

  public void validate(Symbol s)throws IllegalSymbolException{
    if(! contains(s)){
      throw new IllegalSymbolException(
          s, s.getName()+" is not a valid part of "+getName());
    }
  }

  /**
   * Getter for the <code>MaskingDetector<code>
   * @return the <code>MaskingDetector<code>
   */
  public MaskingDetector getMaskingDetector(){
    return maskingDetector;
  }

  public SymbolTokenization getTokenization(String type)
      throws BioException{
    return new CaseSensitiveTokenization(this, type);
  }

  public int size(){
      return delegateAlpha.size();
  }

  public Iterator iterator(){
    return delegateAlpha.iterator();
  }

  /**
   * <code>SoftMaskedAlphabet</code>s cannot add new <code>Symbol</code>s. A
   * <code>ChangeVetoException</code> will be thrown.
   * @param s the <code>Symbol</code> to add.
   * @throws ChangeVetoException when called.
   */
  public void addSymbol(Symbol s) throws ChangeVetoException{
    throw new ChangeVetoException("SoftMaskedAlphabets cannot add new Symbols");
  }

  /**
   * <code>SoftMaskedAlphabet</code>s cannot remove <code>Symbol</code>s. A
   * <code>ChangeVetoException</code> will be thrown.
   * @param s the <code>Symbol</code> to remove.
   * @throws ChangeVetoException when called.
   */
  public void removeSymbol(Symbol s) throws ChangeVetoException{
    throw new ChangeVetoException("SoftMaskedAlphabets cannot remove Symbols");
  }

  /**
   * Determines if a <code>Symbol</code> is masked.
   * @return true if <code>s</code> is masked.
   * @param s the <code>Symbol</code> to test.
   */
  public boolean isMasked (BasisSymbol s) throws IllegalSymbolException {
    validate(s);

    IntegerSymbol b = (IntegerSymbol)s.getSymbols().get(1);
    return (b.intValue() == 1);
  }

  /**
   * Implementations will define how soft masking looks. The
   * <code>DEFAULT</code> implementation considers softmasking to be represented
   * by lower case characters.
   *
   * <p>Copyright (c) 2004 Novartis Institute for Tropical Diseases</p>
   * @author Mark Schreiber
   * @version 1.0
   */
  public interface MaskingDetector{
    public boolean isMasked (String token);

    /**
     * Present the token for a <code>Symbol</code> as it would appear if masked
     * @param token the <code>String</code> to mask.
     * @return the masked token
     */
    public String mask (String token);

    /**
     * Present the token for a <code>Symbol</code> as it would appear if
     * it wasn't softmasked
     * @param token the <code>String</code> to un-mask.
     * @return the un-masked token
     */
    public String unmask (String token);
    public static MaskingDetector DEFAULT = new DefaultMaskingDetector();

    class DefaultMaskingDetector implements MaskingDetector{

      /**
       * Default Behaivour is that if the whole token is lower case it is
       * masked.
       * @param token the <code>String</code> to check for masking
       * @return true is it is all lower case, otherwise false.
       */
      public boolean isMasked(String token){

        for (int i = 0; i < token.length(); i++) {
          if(Character.isUpperCase(token.charAt(i))){
            return false;
          }
        }

        return true;
      }

      /**
       * Masks a token by making it lowercase
       * @param token the <code>String</code> to mask
       * @return a lower case <code>String</code>
       */
      public String mask(String token){
        return token.toLowerCase();
      }

      /**
       * Un-masks the token by making it upper case.
       * @param token the <code>String</code> to unmask
       * @return the upper case <code>String</code>
       */
      public String unmask(String token){
        return token.toUpperCase();
      }
    }
  }

  /**
   * This <code>SymbolTokenizer</code> works with a delegate to softmask
   * symbol tokenization as appropriate. It should only be used in combination
   * with a SoftMaskedAlphabet.
   * You will never instantiate one of these yourself.
   *
   * <p> Copyright (c) 2004 Novartis Institute for Tropical Diseases</p>
   * @author Mark Schreiber
   * @version 1.0
   */
  public class CaseSensitiveTokenization
      extends Unchangeable implements SymbolTokenization{

    private SymbolTokenization delegate;
    private SoftMaskedAlphabet alpha;

    private CaseSensitiveTokenization(
        SoftMaskedAlphabet alpha, String type)
        throws BioException{

      this.alpha = alpha;
      this.delegate = alpha.getMaskedAlphabet().getTokenization(type);
    }

    public Annotation getAnnotation(){
      return Annotation.EMPTY_ANNOTATION;
    }

    public Alphabet getAlphabet(){
      return alpha;
    }

    public SymbolTokenization.TokenType getTokenType(){
      return delegate.getTokenType();
    }

    public Symbol parseToken(String token) throws IllegalSymbolException{
      MaskingDetector md = alpha.getMaskingDetector();
      IntegerSymbol bin;

      Symbol component = delegate.parseToken(token);

      if(md.isMasked(token)){
        bin = binary.getSymbol(1);
      }else{
        bin = binary.getSymbol(0);
      }

      return alpha.getSymbol(new ListTools.Doublet(component, bin));
    }

    public String tokenizeSymbolList(SymbolList sl) throws
        IllegalSymbolException {

      StringBuffer sb = new StringBuffer(sl.length());
      for(int i = 1; i <= sl.length(); i++){
        sb.append(tokenizeSymbol(sl.symbolAt(i)));
      }
      return sb.toString();
    }

    /**
     * The current implementation only supports character parsing. Word or
     * fixed width parsing is not yet supported.
     *
     * @param l the <code>SeqIOListener</code> to callback to.
     * @return a <code>StreamParser</code> that the <code>SeqIOListener</code>
     * talks to.
     */
    public StreamParser parseStream(SeqIOListener l){
      return new CharStreamParser(l);
    }

    public String tokenizeSymbol (Symbol s) throws IllegalSymbolException{
      validate(s);
      Symbol a = (Symbol) ((BasisSymbol)s).getSymbols().get(0);
      String token = delegate.tokenizeSymbol(a);

      if(alpha.isMasked((BasisSymbol) s)){
        return maskingDetector.mask(token);
      }

      return maskingDetector.unmask(token);
    }

    private class CharStreamParser implements StreamParser {
        private SeqIOListener listener;
        private Symbol[] buffer;

        public CharStreamParser(SeqIOListener l) {
            this.listener = l;
            buffer = new Symbol[256];
        }

        public void characters(char[] data, int start, int len)
            throws IllegalSymbolException{
            int cnt = 0;
            while (cnt < len) {
                int bcnt = 0;
                while (cnt < len && bcnt < buffer.length) {
                    buffer[bcnt++] = parseToken(
                      new String(""+data[start + (cnt++)]));
                }
                try {
                    listener.addSymbols(getAlphabet(),
                                        buffer,
                                        0,
                                        bcnt);
                } catch (IllegalAlphabetException ex) {
                    throw new BioError(
                      "Assertion failed: can't add symbols.", ex);
                }
            }
        }

        public void close() {
        }
    }

  }
}
