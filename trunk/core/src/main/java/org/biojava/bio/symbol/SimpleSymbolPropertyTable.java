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
 * Created on December 20, 2000, 7:19 PM
 */

package org.biojava.bio.symbol;

import java.util.HashMap;
import java.util.Map;

/**
 * Class that implements the SymbolPropertyTable interface
 * @author  Mike Jones (primary author)
 * @author  David Huen (minor, extended class to cover pK)
 */
public class SimpleSymbolPropertyTable implements SymbolPropertyTable {

    //Finite ? 
    private final Alphabet source;
    private String name;
    
    private final Map doublePropMap;
    
    public SimpleSymbolPropertyTable(Alphabet source, String name)
    {
        this.source = source;
        this.name = name;
        doublePropMap = new HashMap();
    }
    
    public void setDoubleProperty(Symbol s, String value)
            throws IllegalSymbolException {
                        
        source.validate(s);
        doublePropMap.put(s, new Double(value));
    }

    public String getName() {
        return name;
    }
    
    public Alphabet getAlphabet() {
        return source;
    }
  
    public double getDoubleValue(Symbol s) throws IllegalSymbolException {
        source.validate(s);
        
        Double  value = (Double)doublePropMap.get(s);
        if(value==null) {
            throw new IllegalSymbolException(
                    "Symbol " + s.getName() + " not found in property table " + this.getName());
        }
        return  value.doubleValue();
    }
}
