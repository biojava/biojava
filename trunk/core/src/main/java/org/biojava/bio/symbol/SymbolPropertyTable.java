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
 * Created on December 20, 2000, 7:15 PM
 */

package org.biojava.bio.symbol;

/**
 * class for maintaining properties associated with a symbol
 * @author Mike Jones
 * @author George Waldon
 */
public interface SymbolPropertyTable {
  
  //amino acid mass properties
  public static String AVG_MASS = "avgMass";

  public static String MONO_MASS = "monoMass";

  //amino acid pK properties
  public static String PK_Nterm = "pK_Nterm";
  
  public static String PK = "pK";
  
  public static String PK_Cterm = "pK_Cterm";

  //amino acid Hydropathicity properties
  public static String HYDROPATHICITY = "hydropathicity";
  
  // the name of the property e.g. "isotopic mass"
  public String getName();

  // the alphabet that this property is defined for e.g. PROTEIN
  public Alphabet getAlphabet();

  // the value of the property for a given symbol
  public double getDoubleValue(Symbol s) throws IllegalSymbolException;

 // public void setDoubleProperty(Symbol s, String value) throws IllegalSymbolException;

  // the value of the property for a given symbol
 // public void setDoubleProperty(Symbol s, String value) throws IllegalSymbolException;

}
