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
package org.biojava3.core.symbol;

import java.io.Serializable;

/**
 * A tagged symbol is a simple pair of symbol and tag objects.
 * @author Richard Holland
 * @since 3.0
 * @param T the type of object for the tag.
 */
public interface TaggedSymbol<T> extends Serializable {

    /**
     * Get the symbol for this pair.
     * @return the symbol.
     */
    public Symbol getSymbol();

    /**
     * Get the tag for this pair.
     * @return the tag.
     */
    public T getTag();

    /**
     * Set the symbol for this pair.
     * @param sym the symbol.
     * @return the symbol that used to be in this pair.
     */
    public Symbol setSymbol(Symbol sym);

    /**
     * Set the tag for this pair.
     * @param tah the symbol.
     * @return the tag that used to be in this pair.
     */
    public T setTag(T tag);
}
