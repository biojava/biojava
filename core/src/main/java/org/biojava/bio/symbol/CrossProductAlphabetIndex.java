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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeVetoException;

/**
 * Implementation of AlphabetIndex that handles CrossProductAlphabets
 *
 * @author Kalle N�slund
 * @author Matthew Pocock
 * @author David Huen
 *
 * @since 1.3
 */

// things to change and improve perhaps ?
// 1) need to unregister changelistener at destruction otherwise the Alphabet
//    we index might hold reference and this object will never be garbage
//    collected ?
// 2) just dont rethrow the IndexOutOfBounds exception in the
//    symbolForIndex method

class CrossProductAlphabetIndex extends AbstractChangeable
        implements AlphabetIndex, ChangeListener, Serializable {
    // The alphabet we are indexing
    final FiniteAlphabet	Alpha;
    // list holding AlphabetIndexes for the alphabets making upp the alphabet we index
    List		alphaIndexes;
    List    revIndexes;
    
    public CrossProductAlphabetIndex( FiniteAlphabet theAlpha ) {
        // the alpha we index over
        Alpha = theAlpha;
        // listen to changes in our alpha
        Alpha.addChangeListener( this, Alphabet.SYMBOLS );
        // we need this to compute the index
        alphaIndexes = buildIndexList( theAlpha );
        revIndexes = new ArrayList(alphaIndexes);
        Collections.reverse(revIndexes);
    }
    
    public FiniteAlphabet getAlphabet() {
        return Alpha;
    }
    
    /**
     * This just takes the alphabet, and splits it up into its sub alphabets
     * then gets an AlphabetIndex for each sub alphabet and stores that in an
     * List. This list is then used for the index computation
     */
    private List buildIndexList( FiniteAlphabet fA ) {
        List	subAlphas 	= fA.getAlphabets();
        List	retList		= new ArrayList();
        
        for( Iterator it = subAlphas.iterator() ; it.hasNext() ; ) {
            // no need to check its FiniteAlphabets as we start with a FiniteAlphabet
            FiniteAlphabet currAlpha = ( FiniteAlphabet )it.next();
            AlphabetIndex currInd = AlphabetManager.getAlphabetIndex( currAlpha );
            retList.add( currInd );
        }
        //System.out.println( retList );
        return retList;
    }
    
    public void preChange( ChangeEvent cE ) throws ChangeVetoException {
        // Should perhaps ask eventual listeners to this object if they think
        // the change is ok, and refuse if they dont think it is ????????????????
    }
    
    public void postChange( ChangeEvent cE ) {
        // ok, alphabet have changed, so we better rebuild or list of AlphabetIndexes
        alphaIndexes = buildIndexList( Alpha );
        revIndexes = new ArrayList(alphaIndexes);
        Collections.reverse(revIndexes);
    }
    
    public int indexForSymbol( Symbol s ) throws IllegalSymbolException {
        // ok, split the symbol up into its factors
        List factors = ( ( BasisSymbol )s ).getSymbols();
        if(factors.size() != getAlphabet().getAlphabets().size()) {
            getAlphabet().validate(s);
        }
        
        int	index	= 0; // this is the index we want to compute
        
        // iterate over each factor symbol, and at the same time iterate over the AlphabetIndexList
        // and feed each symbol to its corresponding AlphabetIndex this should perhasp be iteration
        // via arrays instead ? if that speeds it up
        Iterator indexIt  = alphaIndexes.iterator();
        Iterator symbolIt = factors.iterator();
        
        try {
            while( symbolIt.hasNext() ) {
                Symbol 		currentSymbol 	= ( Symbol )symbolIt.next();
                AlphabetIndex 	currentAlphaInd = ( AlphabetIndex )indexIt.next();
                FiniteAlphabet	currentAlphabet	= currentAlphaInd.getAlphabet();
                
                index = index * currentAlphabet.size() + currentAlphaInd.indexForSymbol( currentSymbol );
            }
        } catch (IllegalSymbolException ise) {
            getAlphabet().validate(s);
        }
        
        if(symbolIt.hasNext()) {
            getAlphabet().validate(s);
            throw new BioError("Assertion failure: Ran out of indexers for symbols");
        }
        
        return index;
    }
    
    public Symbol symbolForIndex( int index ) throws IndexOutOfBoundsException {
        List	symbols	= new ArrayList();
        try {
            for ( Iterator it = revIndexes.iterator() ; it.hasNext() ; ) {
                AlphabetIndex	curAlpIndex 	= ( AlphabetIndex ) it.next();
                int		curAlpSize 	= curAlpIndex.getAlphabet().size();
                int		curSymIndValue 	= index % curAlpSize;
                index = index / curAlpSize;
                symbols.add( 0, curAlpIndex.symbolForIndex( curSymIndValue ) );
            }
            
            return( Alpha.getSymbol( symbols ) );
        } catch( IllegalSymbolException isE ) {
            throw new BioError(isE);
        } catch( IndexOutOfBoundsException ioobE ) {
            // this is most likely bad, but i just rethrow the exception, i should do something
            // clever here
            throw ioobE;
        }
        
    }
}
