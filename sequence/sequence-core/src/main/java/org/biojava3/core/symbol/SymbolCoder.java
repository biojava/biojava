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
 * The purpose of the symbol coder is to convert symbols to/from a textual 
 * representation.
 * @author Richard Holland
 * @since 3.0
 */
public interface SymbolCoder extends Serializable {

    /**
     * Decode a string into a symbol. The method must NOT return {@code null}
     * or throw any exception other than {@link NullPointerException} (if the
     * string passed in is null).
     * @param str the string.
     * @return the symbol.
     */
    public Symbol decodeSymbol(String str);

    /**
     * Encode the symbol into a string. The method must NOT return {@code null}
     * or throw any exception other than {@link NullPointerException} (if the
     * symbol passed in is null).
     * @param sym the symbol.
     * @return the string.
     */
    public String encodeSymbol(Symbol sym);

    /**
     * All implementations should extend this implementation and delegate back
     * to it if they are unable to peform the conversion for any reason at all. 
     */
    public static class DefaultSymbolCoder implements SymbolCoder {
		
		private static final long serialVersionUID = 1L;

        public String encodeSymbol(Symbol sym) {
            return sym.toString();
        }

        public Symbol decodeSymbol(String str) {
            return Symbol.get(str);
        }
    }

    /**
     * Used for converting symbols to/from strings.
     */
    public static class StringSymbolCoder extends DefaultSymbolCoder {		
		private static final long serialVersionUID = 1L;
   }

    /**
     * Used for converting symbols to/from characters.
     */
    public static class CharacterSymbolCoder extends DefaultSymbolCoder {
		
		private static final long serialVersionUID = 1L;

        @Override
        public Symbol decodeSymbol(String str) {
            if (str.length() == 1) {
                return Symbol.get(str.charAt(0));
            }
            return super.decodeSymbol(str);
        }
    }

    /**
     * Used for converting symbols to/from integers.
     */
    public static class IntegerSymbolCoder extends DefaultSymbolCoder {
		
		private static final long serialVersionUID = 1L;

        @Override
        public Symbol decodeSymbol(String str) {
            try {
                return Symbol.get(Integer.valueOf(str));
            } catch (NumberFormatException nfe) {
            // Ignore and fall back to default.
            }
            return super.decodeSymbol(str);
        }
    }

    /**
     * Used for converting symbols to/from doubles.
     */
    public static class DoubleSymbolCoder extends DefaultSymbolCoder {
		
		private static final long serialVersionUID = 1L;

        @Override
        public Symbol decodeSymbol(String str) {
            try {
                return Symbol.get(Double.valueOf(str));
            } catch (NumberFormatException nfe) {
            // Ignore and fall back to default.
            }
            return super.decodeSymbol(str);
        }
    }

    /**
     * Used for converting symbols to/from booleans.
     */
    public static class BooleanSymbolCoder extends DefaultSymbolCoder {
		
		private static final long serialVersionUID = 1L;

        @Override
        public Symbol decodeSymbol(String str) {
            if (str.length() == 1) {
                if ("0".equals(str)) {
                    return Symbol.get(Boolean.FALSE);
                }
                if ("1".equals(str)) {
                    return Symbol.get(Boolean.TRUE);
                }
            }
            return super.decodeSymbol(str);
        }

        @Override
        public String encodeSymbol(Symbol sym) {
            if (sym.getObject() instanceof Boolean) {
                return (((Boolean) sym.getObject()).booleanValue()) ? "1" : "0";
            }
            return super.encodeSymbol(sym);
        }
    }
}
