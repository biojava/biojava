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

import junit.framework.TestCase;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ListTools;


public class TestSoftMaskedAlphabet extends TestCase {
  private SoftMaskedAlphabet softMaskedAlphabet = null;
  private AtomicSymbol masked;
  private AtomicSymbol unmasked;

  protected void setUp() throws Exception {
    super.setUp();
    softMaskedAlphabet = SoftMaskedAlphabet.getInstance(DNATools.getDNA());
    masked = (AtomicSymbol)softMaskedAlphabet.getSymbol(
        new ListTools.Doublet(DNATools.a(),
                             IntegerAlphabet.getInstance().getSymbol(1)));
   unmasked = (AtomicSymbol)softMaskedAlphabet.getSymbol(
       new ListTools.Doublet(DNATools.a(),
                            IntegerAlphabet.getInstance().getSymbol(0)));
  }

  protected void tearDown() throws Exception {
    softMaskedAlphabet = null;
    masked = null;
    unmasked = null;
    super.tearDown();
  }

  public void testAddSymbol(){
    Symbol s = DNATools.a();
    try{
      softMaskedAlphabet.addSymbol(s);
    }catch(ChangeVetoException ex){
      return;
    }
    fail("Should have thrown an exception");
  }

  public void testContains() {
    assertTrue(softMaskedAlphabet.contains(masked));
    assertTrue(softMaskedAlphabet.contains(unmasked));
  }

  public void testGetAlphabets() {
    Alphabet expectedReturn = null;
    Alphabet actualReturn = null;

    List l = softMaskedAlphabet.getAlphabets();
    assertTrue(l.size() == 2);

    expectedReturn = DNATools.getDNA();
    assertTrue(l.get(0) instanceof Alphabet);
    actualReturn = (Alphabet)l.get(0);
    assertEquals("return value", expectedReturn, actualReturn);

    expectedReturn = IntegerAlphabet.getSubAlphabet(0,1);
    assertTrue(l.get(1) instanceof Alphabet);
    actualReturn = (Alphabet)l.get(1);
    assertEquals("return value", expectedReturn, actualReturn);

  }

  public void testGetAmbiguity(){
    Set s = null;
    try{
      Symbol actualReturn = softMaskedAlphabet.getAmbiguity(s);
      actualReturn=actualReturn==null?null:actualReturn; // trick
    }catch(UnsupportedOperationException ex){
      return;
    }
    fail("Should have thrown and UnsupportedOperationException");
  }

  public void testGetAnnotation() {
    assertNotNull(softMaskedAlphabet.getAnnotation());
  }

  public void testGetGapSymbol() {
    Symbol expectedReturn = AlphabetManager.getGapSymbol(
        softMaskedAlphabet.getAlphabets());
    Symbol actualReturn = softMaskedAlphabet.getGapSymbol();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetMaskedAlphabet() {
    FiniteAlphabet expectedReturn = DNATools.getDNA();
    FiniteAlphabet actualReturn = softMaskedAlphabet.getMaskedAlphabet();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetMaskingDetector() {
    SoftMaskedAlphabet.MaskingDetector expectedReturn =
        SoftMaskedAlphabet.MaskingDetector.DEFAULT;
    SoftMaskedAlphabet.MaskingDetector actualReturn =
        softMaskedAlphabet.getMaskingDetector();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetName() {
    String expectedReturn = "Softmasked {DNA}";
    String actualReturn = softMaskedAlphabet.getName();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetSymbol() throws IllegalSymbolException {
    List l = new ListTools.Doublet(DNATools.g(),
                           IntegerAlphabet.getInstance().getSymbol(1));
    Symbol expectedReturn =
        softMaskedAlphabet.getDelegate().getSymbol(l);
    Symbol actualReturn = softMaskedAlphabet.getSymbol(l);
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testGetTokenization() throws BioException {
    String type = "token";
    SymbolTokenization actualReturn = softMaskedAlphabet.getTokenization(type);
    assertNotNull(actualReturn);
    assertTrue(actualReturn instanceof SoftMaskedAlphabet.CaseSensitiveTokenization);
  }

  public void testTokenizeSymbol() throws BioException{
    SymbolTokenization st = softMaskedAlphabet.getTokenization("token");
    assertEquals(softMaskedAlphabet, st.getAlphabet());

    AtomicSymbol s = null;
    s = (AtomicSymbol)st.parseToken("a");
    assertTrue(softMaskedAlphabet.isMasked(s));
    assertEquals("a", st.tokenizeSymbol(s));

    s = (AtomicSymbol) st.parseToken("c");
    assertTrue(softMaskedAlphabet.isMasked(s));
    assertEquals("c", st.tokenizeSymbol(s));

    s = (AtomicSymbol) st.parseToken("g");
    assertTrue(softMaskedAlphabet.isMasked(s));
    assertEquals("g", st.tokenizeSymbol(s));

    s = (AtomicSymbol) st.parseToken("t");
    assertTrue(softMaskedAlphabet.isMasked(s));
    assertEquals("t", st.tokenizeSymbol(s));

    s = (AtomicSymbol)st.parseToken("A");
    assertTrue(! softMaskedAlphabet.isMasked(s));
    assertEquals("A", st.tokenizeSymbol(s));

    s = (AtomicSymbol) st.parseToken("C");
    assertTrue(! softMaskedAlphabet.isMasked(s));
    assertEquals("C", st.tokenizeSymbol(s));

    s = (AtomicSymbol) st.parseToken("G");
    assertTrue(! softMaskedAlphabet.isMasked(s));
    assertEquals("G", st.tokenizeSymbol(s));

    s = (AtomicSymbol) st.parseToken("T");
    assertTrue(! softMaskedAlphabet.isMasked(s));
    assertEquals("T", st.tokenizeSymbol(s));
  }

  public void testTokenization() throws IllegalSymbolException,
      BioException {

    String seqString = "actgATGC";
    SymbolList sl = new SimpleSymbolList(
        softMaskedAlphabet.getTokenization("token"),
        seqString);

    assertEquals(seqString, sl.seqString());
  }

  public void testIsMasked() throws IllegalSymbolException {
    List a = new ListTools.Doublet(DNATools.g(),
                           IntegerAlphabet.getInstance().getSymbol(1));
    List b = new ListTools.Doublet(DNATools.g(),
                           IntegerAlphabet.getInstance().getSymbol(0));

    AtomicSymbol s = (AtomicSymbol) softMaskedAlphabet.getSymbol(a);
    boolean expectedReturn = true;
    boolean actualReturn = softMaskedAlphabet.isMasked(s);
    assertEquals("return value", expectedReturn, actualReturn);

    s = (AtomicSymbol) softMaskedAlphabet.getSymbol(b);
    expectedReturn = false;
    actualReturn = softMaskedAlphabet.isMasked(s);
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testIterator() {
    Iterator i = softMaskedAlphabet.iterator();
    int count = 0;
    while(i.hasNext()){
      assertTrue(i.next() instanceof AtomicSymbol);
      count ++;
    }
    assertTrue(count == 8);
  }

  public void testRemoveSymbol(){
    Symbol s = null;
    try {
      softMaskedAlphabet.removeSymbol(s);
    }
    catch (ChangeVetoException ex) {
      return;
    }
    fail("Should have thrown ChangeVetoException");
  }

  public void testSize() {
    int expectedReturn = 8;
    int actualReturn = softMaskedAlphabet.size();
    assertEquals("return value", expectedReturn, actualReturn);
  }

  public void testValidate() {
    for (Iterator iter = softMaskedAlphabet.iterator(); iter.hasNext(); ) {
      Symbol s = (Symbol) iter.next();
      try {
        softMaskedAlphabet.validate(s);
      }
      catch (IllegalSymbolException ex) {
        fail("Symbol "+s.getName()+" not a valid member of this alphabet");
      }
    }
  }
}
