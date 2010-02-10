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

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.ManyToOneTranslationTable;
import org.biojava.bio.symbol.SimpleGeneticCodeTable;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.TranslationTable;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <code>RNAToolsTest</code> tests are to ensure that the static methods
 * provided by <code>RNATools</code> work as advertised (including under
 * error conditions) and are internally consistent.
 *
 * The classes to which <code>RNATools</code> delegates e.g.
 * <code>AlphabetManager</code> and <code>SymbolTokenization</code> should
 * be tested in their own unit tests.
 *
 * @author Moses Hohman
 * @author gwaldon
 */
public class RNAToolsTest extends TestCase {
    public RNAToolsTest(String name){
      super(name);
    }
    public void testSymbols() throws IllegalSymbolException {
        assertEquals("a", getRNATokenization().tokenizeSymbol(RNATools.a()));
        assertEquals("c", getRNATokenization().tokenizeSymbol(RNATools.c()));
        assertEquals("g", getRNATokenization().tokenizeSymbol(RNATools.g()));
        assertEquals("u", getRNATokenization().tokenizeSymbol(RNATools.u()));
    }

    public void testCreateRNASequenceGetName() throws IllegalSymbolException {
        String sequenceName = "the name";
        Sequence sequence = RNATools.createRNASequence("uacg-", sequenceName);
        assertEquals(sequenceName, sequence.getName());
    }

    public void testCreateRNASequenceAnnotationIsEmpty() throws IllegalSymbolException {
        Sequence sequence = RNATools.createRNASequence("uacg-", "the name");
        assertTrue(sequence.getAnnotation().asMap().isEmpty());
    }

    public void testCreateRNASequenceURNIsBlank() throws IllegalSymbolException {
        Sequence sequence = RNATools.createRNASequence("uacg-", "the name");
        assertEquals("", sequence.getURN());
    }

    public void testCreateRNASequenceSequenceCorrect() throws IllegalSymbolException {
        String sequenceString = "uacg-";
        assertEquals(sequenceString, ((SimpleSequence) RNATools.createRNASequence(sequenceString, "the name")).seqString());
    }

    public void testGetCodonAlphabet() {
        List threeCopiesOfRNAAlphabet = Collections.nCopies(3, RNATools.getRNA());
        assertEquals(RNATools.getCodonAlphabet(), AlphabetManager.getCrossProductAlphabet(threeCopiesOfRNAAlphabet));
    }

    public void testCreateRNA() throws IllegalAlphabetException, IllegalSymbolException {
        assertEquals("uacg-", getRNATokenization().tokenizeSymbolList(RNATools.createRNA("uacg-")));
    }

    public void testComplementA() throws IllegalSymbolException {
        assertEquals(RNATools.u(), RNATools.complement(RNATools.a()));
    }

    public void testComplementC() throws IllegalSymbolException {
        assertEquals(RNATools.g(), RNATools.complement(RNATools.c()));
    }

    public void testComplementG() throws IllegalSymbolException {
        assertEquals(RNATools.c(), RNATools.complement(RNATools.g()));
    }

    public void testComplementU() throws IllegalSymbolException {
        assertEquals(RNATools.a(), RNATools.complement(RNATools.u()));
    }

    public void testComplementSymbolFromAnotherAlphabet() {
        try {
            RNATools.complement(DNATools.t());
            fail("Should have thrown an IllegalSymbolException");
        } catch (IllegalSymbolException e) {
        }
    }

    public void testComplementBadSymbol() throws IllegalSymbolException {
        try {
            RNATools.complement(createBadSymbol());
            fail("Should have thrown a BioError");
        } catch (BioError e) {
        }
    }

    public void testComplementList() throws IllegalAlphabetException, IllegalSymbolException {
        assertEquals(RNATools.createRNA("augc"), RNATools.complement(RNATools.createRNA("uacg")));
    }

    public void testReverseComplement() throws IllegalAlphabetException, IllegalSymbolException {
        assertEquals(RNATools.createRNA("cgua"), RNATools.reverseComplement(RNATools.createRNA("uacg")));
    }

    public void testTranscribe() throws IllegalAlphabetException, IllegalSymbolException {
        assertEquals(RNATools.createRNA("gcau"), RNATools.transcribe(DNATools.createDNA("gcat")));
    }

    public void testTranslate() throws IllegalAlphabetException, IllegalSymbolException {
        assertEquals(ProteinTools.createProtein(
                "FFLLSSSSYY**CC*W"
                + "LLLLPPPPHHQQRRRR"
                + "IIIMTTTTNNKKSSRR"
                + "VVVVAAAADDEEGGGG"),
                RNATools.translate(RNATools.createRNA(
                        "uuuuucuuauug" + "ucuuccucaucg" + "uauuacuaauag" + "uguugcugaugg"
                + "cuucuccuacug" + "ccucccccaccg" + "caucaccaacag" + "cgucgccgacgg"
                + "auuaucauaaug" + "acuaccacaacg" + "aauaacaaaaag" + "aguagcagaagg"
                + "guugucguagug" + "gcugccgcagcg" + "gaugacgaagag" + "gguggcggaggg"
                )));
    }

    public void testForSymbolA() throws IllegalSymbolException {
        assertEquals(RNATools.a(), RNATools.forSymbol('a'));
    }

    public void testForSymbolC() throws IllegalSymbolException {
        assertEquals(RNATools.c(), RNATools.forSymbol('c'));
    }

    public void testForSymbolG() throws IllegalSymbolException {
        assertEquals(RNATools.g(), RNATools.forSymbol('g'));
    }

    public void testForSymbolU() throws IllegalSymbolException {
        assertEquals(RNATools.u(), RNATools.forSymbol('u'));
    }

    public void testForSymbolWithBadChar() {
        try {
            RNATools.forSymbol('$');
            fail("Should have thrown an IllegalSymbolException");
        } catch (IllegalSymbolException e) {
        }
    }

    public void testForIndex0() throws IllegalSymbolException {
        assertEquals(RNATools.a(), RNATools.forIndex(0));
    }

    public void testForIndex1() throws IllegalSymbolException {
        assertEquals(RNATools.g(), RNATools.forIndex(1));
    }

    public void testForIndex2() throws IllegalSymbolException {
        assertEquals(RNATools.c(), RNATools.forIndex(2));
    }

    public void testForIndex3() throws IllegalSymbolException {
        assertEquals(RNATools.u(), RNATools.forIndex(3));
    }

    public void testForIndexWithBadIndex() {
        try {
            RNATools.forIndex(4);
            fail("Should have thrown an IndexOutOfBoundsException");
        } catch (IndexOutOfBoundsException e) {
        }
    }

    public void testIndexA() throws IllegalSymbolException {
        assertEquals(0, RNATools.index(RNATools.a()));
    }

    public void testIndexG() throws IllegalSymbolException {
        assertEquals(1, RNATools.index(RNATools.g()));
    }

    public void testIndexC() throws IllegalSymbolException {
        assertEquals(2, RNATools.index(RNATools.c()));
    }

    public void testIndexU() throws IllegalSymbolException {
        assertEquals(3, RNATools.index(RNATools.u()));
    }

    public void testIndexSymbolFromAnotherAlphabet() {
        try {
            RNATools.index(DNATools.t());
            fail("Should have thrown an IllegalSymbolException");
        } catch (IllegalSymbolException e) {
        }
    }

    public void testIndexBadSymbol() {
        try {
            RNATools.index(createBadSymbol());
            fail("Should have thrown an IllegalSymbolException");
        } catch (IllegalSymbolException e) {
        }
    }

    public void testGeneticCode() {
        try {
            // get the universal genetic code and give it a trashing
            ManyToOneTranslationTable geneticCode = RNATools.getGeneticCode(TranslationTable.UNIVERSAL);
            assertNotNull(geneticCode);

            // for the entire amino-acid alphabet, do the reverse lookup
            // and then confirm that the Set of Symbols return is correct
            FiniteAlphabet aaAlfa = ProteinTools.getTAlphabet();
            assertNotNull(aaAlfa);

            for (Iterator aaI = aaAlfa.iterator(); aaI.hasNext();) {
                Symbol residue = (Symbol) aaI.next();

                // get the List of codons that yield this amino-acid
                Set codons = geneticCode.untranslate(residue);
                assertNotNull(codons);

                // iterate thru the list confirming that they correspond
                // to the expected amino-acid
                for (Iterator codonI = codons.iterator(); codonI.hasNext(); ) {
                    Symbol codon = (Symbol) codonI.next();

                    // translate codon
                    Symbol xlatedResidue = geneticCode.translate(codon);

                    assertEquals(residue, xlatedResidue);
                }
            }
            
            // The following test is dependant of the DDBJ/EMBL/GenBank Feature Table,
            // Version 6.5  Apr 2006.
            // get all table names
            Set tableNames = RNATools.getGeneticCodeNames();
            assertEquals(17,tableNames.size()); //total 17 tables
            
            //test standard table
            ManyToOneTranslationTable geneticCode1 = RNATools.getGeneticCode(TranslationTable.UNIVERSAL);
            ManyToOneTranslationTable geneticCode2 = RNATools.getGeneticCode(1);
            assertEquals((Object)geneticCode1, (Object)geneticCode2);
            String description = ((SimpleGeneticCodeTable)geneticCode2).getDescription();
            assertEquals(description,"Standard Code");
            int number = ((SimpleGeneticCodeTable)geneticCode2).getTableNumber();
            assertEquals(1,number);
     
        }
        catch (IllegalSymbolException ise) {
        }
    }

    public void assertEquals(Symbol expected, Symbol actual) throws IllegalSymbolException {
        assertEquals(getRNATokenization().tokenizeSymbol(expected), getRNATokenization().tokenizeSymbol(actual));
    }

    public void assertEquals(SymbolList expected, SymbolList actual) throws IllegalAlphabetException, IllegalSymbolException {
        try {
            SymbolTokenization expectedTokenization = expected.getAlphabet().getTokenization("token");
            SymbolTokenization actualTokenization = actual.getAlphabet().getTokenization("token");
            assertEquals(expectedTokenization.tokenizeSymbolList(expected), actualTokenization.tokenizeSymbolList(actual));
        } catch (BioException e) {
            throw new BioError(e);
        }
    }

    private SymbolTokenization getRNATokenization() {
        try {
            return RNATools.getRNA().getTokenization("token");
        } catch (BioException e) {
            throw new BioError(e);
        }
    }

    private Symbol createBadSymbol() {
        return new Symbol() {
            public String getName() {
                return null;
            }

            public Alphabet getMatches() {
                return new FiniteAlphabet() {
                    public int size() {
                        return 0;
                    }

                    public String getName() {
                        return null;
                    }

                    public Annotation getAnnotation() {
                        return null;
                    }

                    public void addChangeListener(ChangeListener cl) {
                    }

                    public Iterator iterator() {
                        return null;
                    }

                    public List getAlphabets() {
                        return null;
                    }

                    public void addChangeListener(ChangeListener cl, ChangeType ct) {
                    }

                    public void addSymbol(Symbol s)
                            throws IllegalSymbolException, ChangeVetoException {
                    }

                    public Symbol getSymbol(List rl)
                            throws IllegalSymbolException {
                        return null;
                    }

                    public void removeChangeListener(ChangeListener cl) {
                    }

                    public void removeSymbol(Symbol s)
                            throws IllegalSymbolException, ChangeVetoException {
                    }

                    public Symbol getAmbiguity(Set syms)
                            throws IllegalSymbolException {
                        return null;
                    }

                    public void removeChangeListener(ChangeListener cl, ChangeType ct) {
                    }


                    public Symbol getGapSymbol() {
                        return null;
                    }

                    public boolean isUnchanging(ChangeType ct) {
                        return false;
                    }

                    public boolean contains(Symbol s) {
                        return false;
                    }

                    public void validate(Symbol s) throws IllegalSymbolException {
                    }

                    public SymbolTokenization getTokenization(String name) throws BioException {
                        return null;
                    }
                };
            }

            public Annotation getAnnotation() {
                return null;
            }

            public void addChangeListener(ChangeListener cl) {
            }

            public void addChangeListener(ChangeListener cl, ChangeType ct) {
            }

            public void removeChangeListener(ChangeListener cl) {
            }

            public void removeChangeListener(ChangeListener cl, ChangeType ct) {
            }

            public boolean isUnchanging(ChangeType ct) {
                return false;
            }
        };
    }
}
