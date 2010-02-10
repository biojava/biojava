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

package org.biojavax.bio.seq.io;

import org.biojava.bio.seq.io.SequenceBuilder;
import org.biojava.bio.symbol.SymbolListFactory;


/**
 * Generates RichSequenceBuilder objects.
 * @author Mark Schreiber
 * @author Richard Holland
 * @see RichSequenceBuilder
 * @since 1.5
 */
public class SimpleRichSequenceBuilderFactory implements RichSequenceBuilderFactory {
    private SymbolListFactory fact = null;
    private int threshold = 0;
    
    /** Creates a new instance of SimpleRichSequenceBuilderFactory */
    public SimpleRichSequenceBuilderFactory() {
      this(null, 0);
    }

    /** 
     * Creates a new instance of SimpleRichSequenceBuilderFactory
     * @param fact the factory to use when building the <code>SymbolList</code>.
     */
    public SimpleRichSequenceBuilderFactory(SymbolListFactory fact) {
      this(fact, 0);
    }
    
    /** 
     * Creates a new instance of SimpleRichSequenceBuilderFactory that uses
     * a specified factory for <code>SymbolLists</code> longer than a specified
     * length. Before that a <code>SimpleSymbolListFacotry</code> is used.
     * @param fact the factory to use when building the <code>SymbolList</code>.
     * @param threshold the threshold to exceed before using this factory
     */
    public SimpleRichSequenceBuilderFactory(SymbolListFactory fact, int threshold) {
      this.fact = fact;
      this.threshold = threshold;
    }
    /**
     * {@inheritDoc}
     */
    public SequenceBuilder makeSequenceBuilder() {
        if(this.fact == null){
            return new SimpleRichSequenceBuilder(); 
        }else{
            return new SimpleRichSequenceBuilder(this.fact, this.threshold);
        }
    }
    
}
