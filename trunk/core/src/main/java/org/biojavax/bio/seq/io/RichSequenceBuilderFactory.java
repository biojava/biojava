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

import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.symbol.PackedSymbolListFactory;

/**
 * Simple factory for constructing new RichSequenceBuilder objects.
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */
public interface RichSequenceBuilderFactory extends SequenceBuilderFactory {
    
    /**
     * The value that will be used as a threshold for the <code>THRESHOLD</code>
     * builder. Set to 5000.
     */
    public final static int THRESHOLD_VALUE = 5000;
    
    /**
     * Accessor for the default factory. This implementation will not 
     * do any compression of a sequence regardless of size.
     */
    public final static RichSequenceBuilderFactory FACTORY = new SimpleRichSequenceBuilderFactory();
    
    /**
     * Accessor for a factory that produces builders that compress the
     * <code>SymbolList</code> of a <code>RichSequence</code>.
     */
    public final static RichSequenceBuilderFactory PACKED = new SimpleRichSequenceBuilderFactory(new PackedSymbolListFactory(), 0);
    
    /**
     * Accessor for a factory that produces builders that compress the
     * <code>SymbolList</code> of a <code>RichSequence</code> when the length of the
     * <code>SymbolList</code> exceeds <code>THRESHOLD</code>.
     */
    public final static RichSequenceBuilderFactory THRESHOLD = new SimpleRichSequenceBuilderFactory(new PackedSymbolListFactory(), THRESHOLD_VALUE);
    
    
    
}
