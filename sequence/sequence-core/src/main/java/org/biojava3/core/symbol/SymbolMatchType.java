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
 * The various ways in which two symbols can match each other within
 * the confines of a single alphabet.
 * @author Richard Holland
 * @since 3.0
 */
public enum SymbolMatchType implements Serializable {

    /**
     * The symbols do not match.
     */
    MISMATCH,
    /**
     * The symbols match exactly.
     */
    EXACT_MATCH,
    /**
     * The symbols match when converted to strings and compared
     * with {@link String#equals()}.
     */
    EXACT_STRING_MATCH,
    /**
     * The symbols match when converted to strings and compared
     * with {@link String#equalsIgnoreCase()}.
     */
    MIXED_CASE_MATCH,
    /**
     * One of the symbols appears in the other symbol's ambiguity set
     * in this alphabet.
     */
    AMBIGUOUS_MATCH,
    /**
     * One of the symbols appears in the other symbol's ambiguity set
     * in this alphabet but only when compared as a string using 
     * {@link String#equals()},
     * OR one of the symbols appears in the ambiguity set of a symbol
     * in this alphabet that matches the other symbol when compared as 
     * a string using {@link String#equals()},
     */
    AMBIGUOUS_STRING_MATCH,
    /**
     * One of the symbols appears in the other symbol's ambiguity set
     * in this alphabet but only when compared as a string using 
     * {@link String#equalsIgnoreCase()},
     * OR one of the symbols appears in the ambiguity set of a symbol
     * in this alphabet that matches the other symbol when compared as 
     * a string using {@link String#equalsIgnoreCase()},
     */
    AMBIGUOUS_MIXED_CASE_MATCH,
    /**
     * One of the symbols is the gap symbol in this alphabet and the 
     * other is not.
     */
    GAP_MATCH
}
