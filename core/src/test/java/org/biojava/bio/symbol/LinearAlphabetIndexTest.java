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

/*
 * LinearAlphabetIndexTest.java
 * JUnit based test
 *
 * Created on July 10, 2007, 3:19 PM
 */

package org.biojava.bio.symbol;

import junit.framework.*;
import java.lang.ref.Reference;
import java.lang.ref.WeakReference;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 *
 * @author Mark Schreiber
 */
public class LinearAlphabetIndexTest extends TestCase {
    FiniteAlphabet DNA = DNATools.getDNA();
    
    public LinearAlphabetIndexTest(String testName) {
        super(testName);
    }



    public static Test suite() {
        TestSuite suite = new TestSuite(LinearAlphabetIndexTest.class);
        
        return suite;
    }


    /**
     * Test of getAlphabet method, of class org.biojava.bio.symbol.LinearAlphabetIndex.
     */
    public void testGetAlphabet() {
        System.out.println("getAlphabet");
        
        LinearAlphabetIndex instance = new LinearAlphabetIndex(DNA);
        
        FiniteAlphabet expResult = DNA;
        FiniteAlphabet result = instance.getAlphabet();
        assertEquals(expResult, result);
        
    }

    /**
     * Test of symbolForIndex method, of class org.biojava.bio.symbol.LinearAlphabetIndex.
     */
    public void testSymbolForIndex() {
        System.out.println("symbolForIndex");
        
        LinearAlphabetIndex instance = new LinearAlphabetIndex(DNA);
        
        for(int i = 0; i < DNA.size(); i++){
          Symbol result = instance.symbolForIndex(i);
          assertNotNull(result);
          assertTrue(DNA.contains(result));
        }
        //unfortunately we can't test the order of the index as it is not guarenteed
        //to be the same on all JVMs
    }

    /**
     * Test of indexForSymbol method, of class org.biojava.bio.symbol.LinearAlphabetIndex.
     */
    public void testIndexForSymbol() throws Exception {
        System.out.println("indexForSymbol");
        
        for(Iterator<Symbol> i = DNA.iterator(); i.hasNext();){
          Symbol s = i.next();
          LinearAlphabetIndex instance = new LinearAlphabetIndex(DNA);
          int result = instance.indexForSymbol(s);
          assertTrue(result >= 0);
        }
        
    }
    
    /**
     * Tests that the index rebuilds when it's size changes
     */
    public void testBugFix2330() throws Exception{
        SimpleAlphabet alpha = new SimpleAlphabet();
        LinearAlphabetIndex index = new LinearAlphabetIndex(alpha);
        
        Symbol s0 = AlphabetManager.createSymbol("s0");
        alpha.addSymbol(s0);
        assertNotNull(index.symbolForIndex(0));
        assertTrue(index.symbolForIndex(0) == s0);
        assertTrue(index.indexForSymbol(s0) == 0);
        
        
        alpha.removeSymbol(s0);
        try{
            index.symbolForIndex(0);
            fail("Expected exception, there should be no Symbol at index 0");
        }catch (Exception ex){}
        
    }
    
}
