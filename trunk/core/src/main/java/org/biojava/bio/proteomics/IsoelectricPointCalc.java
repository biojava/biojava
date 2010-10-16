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

package org.biojava.bio.proteomics;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolPropertyTable;
import org.biojava.utils.math.BinarySearch;
import org.biojava.utils.math.ComputeObject;

/** Class that computes isoelectric point for denaturated proteins. These pIs 
 * are useful for predicting the position of a protein on a 2D gel.<p>
 * The pK values are taken from Bjellqvist B. et al., "Reference 
 * points for comparisons of two-dimensional maps of proteins from different
 * human cell types defined in a pH scale where isoelectric points correlate
 * with polypeptide compositions", Electrophoresis 1994, 15, 529-539.<p>
 *
 * @author David Huen
 * @author George Waldon
 * @since 1.22
 *
 */
public class IsoelectricPointCalc {
    
   /** minimum pH value */
     private static double PH_MIN = 0.0;
    
    /** maximum pH value */
    private static double PH_MAX = 14.0;
    
    /** desired precision */
    private static double EPSI = 0.0001;

    private static Map pK_NtermCache = new HashMap();
    private static Map pKCache = new HashMap();
    private static Map pK_CtermCache = new HashMap();

    public IsoelectricPointCalc() {
        // recover pK and pK_NTerm tables and cache only relevant residues
        SymbolPropertyTable PK_NtermTable = ProteinTools.getSymbolPropertyTable(SymbolPropertyTable.PK_Nterm);
        SymbolPropertyTable pKTable = ProteinTools.getSymbolPropertyTable(SymbolPropertyTable.PK);
        SymbolPropertyTable PK_CtermTable = ProteinTools.getSymbolPropertyTable(SymbolPropertyTable.PK_Cterm);
        
        Iterator aaSyms = ProteinTools.getAlphabet().iterator();
        
        // iterate thru' all AA symbols and cache the non-zero pKs
        while (aaSyms.hasNext()) {
            Symbol sym = (Symbol) aaSyms.next();
            
            // only cache symbols that have a non-zero pK
            try {
                double pK = PK_NtermTable.getDoubleValue(sym);
                if (Math.abs(pK) > 0.01) {
                    pK_NtermCache.put(sym, new Double(pK));
                }
                
                pK = pKTable.getDoubleValue(sym);
                if (Math.abs(pK) > 0.01) {
                    pKCache.put(sym, new Double(pK));
                }
                
                pK = PK_CtermTable.getDoubleValue(sym);
                if (Math.abs(pK) > 0.01) {
                    pK_CtermCache.put(sym, new Double(pK));
                }
                
            } catch (IllegalSymbolException ise) {
                // SimpleSymbolPropertyTable throws this if there is no value for the symbol
                // just ignore.
            }
        }
    }

    public class ChargeCalculator
            implements ComputeObject {
        Map counts = null;
        Symbol Nterm = null;
        Symbol Cterm = null;
        boolean hasFreeNTerm = true;
        boolean hasFreeCTerm = true;

        private ChargeCalculator(SymbolList peptide, boolean hasFreeNTerm, boolean hasFreeCTerm) {
            counts = residueCount(peptide);
            this.hasFreeNTerm = hasFreeNTerm;
            this.hasFreeCTerm = hasFreeCTerm;
        }

        /**
         * counts up number of times a relevant AA appears in protein.
         */
        private Map residueCount(SymbolList peptide) {
            // iterate thru' peptide collating number of relevant residues
            Iterator residues = peptide.iterator();
            
            Map symbolCounts = new HashMap();
            
            while (residues.hasNext()) {
                Symbol sym = (Symbol) residues.next();
                if(Nterm==null) Nterm = sym;
                if(!residues.hasNext()) Cterm = sym;
                if (pKCache.containsKey(sym)) {
                    // count the residues
                    Integer currCount = (Integer) symbolCounts.get(sym);
                    if (currCount != null) {
                        int currCountAsInt = currCount.intValue();
                        symbolCounts.put(sym, new Integer(++currCountAsInt));
                    } else {
                        symbolCounts.put(sym, new Integer(1));
                    }
                }
            }
            
            return symbolCounts;
        }
        
        /**
         * computes charge at given pH
         */
        public double compute(double pH) {
            double charge = 0.0;
            
            // iterate thru' all counts computing the partial contribution to charge
            Iterator aaI = counts.keySet().iterator();
            
            //by convention negative pK values in ResidueProperties.xml are for
            //acids and inversely for bases.
            while (aaI.hasNext()) {
                // get back the symbol
                Symbol sym = (Symbol) aaI.next();
                
                // retrieve the pK and count
                Double value = (Double) pKCache.get(sym);
                if (value != null) {
                    double pK = value.doubleValue();
                    double count = ((Integer) counts.get(sym)).intValue();
                    boolean isAcid = pK<0;
                    if(isAcid==true) {
                        pK = -pK;
                        double cr = Math.pow(10.0, pK - pH);
                        charge -= count/(cr + 1.0); // -0.5 per aa at pH = pK
                    } else {
                        double cr = Math.pow(10.0, pH - pK);
                        charge += count/(cr + 1.0); // +0.5 per aa at pH = pK
                    }
                }
            }
            
            // N-terminal end charges
            if (hasFreeNTerm) {
                Double value = (Double) pK_NtermCache.get(Nterm);
                if (value != null) {
                    double pK = + value.doubleValue();
                    double cr = Math.pow(10.0, pH - pK);
                    charge += 1/(cr + 1.0);
                }
            }
            
            // C-terminal end charges
            if(hasFreeCTerm) {
                Double value = (Double) pK_CtermCache.get(Cterm);
                if (value != null) {
                    double pK = - value.doubleValue();
                    double cr = Math.pow(10.0, pK - pH);
                    charge -= 1/(cr + 1.0); // -0.5 per aa at pH = pK
                }
            }

            return charge;
        }
    }

    /**
     * Computes isoelectric point of specified peptide.
     *
     * @param peptide peptide of which pI is required.
     * @param hasFreeNTerm has free N-terminal amino group.
     * @param hasFreeCTerm has free C-terminal carboxyl group.
     */
    public double getPI(SymbolList peptide, boolean hasFreeNTerm, boolean hasFreeCTerm)
        throws IllegalAlphabetException, BioException
    {
        // verify that the peptide is really a peptide
        if ( (peptide.getAlphabet() == ProteinTools.getTAlphabet())
        || (peptide.getAlphabet() == ProteinTools.getAlphabet()) ) {
            
            // create object to handle the peptide
            ComputeObject computeObj = new ChargeCalculator(peptide, hasFreeNTerm, hasFreeCTerm);
            
            // solve the charge equation
            double pI = 0.0;
            try {
                pI = BinarySearch.solve(PH_MIN, PH_MAX, EPSI, computeObj);
            } catch( BioException ex) {
                BioException ex2 =  new BioException("Error, the peptide probably contains only positive or negative charges");
                ex2.initCause(ex);
                throw ex2;
            }
            return pI;
        } else {
            // not a peptide
            throw new IllegalAlphabetException();
        }
    }
    
    private static IsoelectricPointCalc calculator;
    
    /** Static public method to compute the pI for a polypeptide in
     * denaturating and reduced conditions with both free ends. 
     * Various ambiguity symbols, symbols for which pK data are not available, or 
     * illegal symbols are not contributing to the calculated pI.<p>
     * This method returns the same values as ExPASy's Compute pI/Mw program
     *
     * @param peptide peptide of which pI is required.
     * @return the calculated pI
     * @since 1.5
     */
    public static double getIsoelectricPoint(SymbolList peptide) 
    throws IllegalAlphabetException, BioException {
        if(calculator==null) {
            calculator = new IsoelectricPointCalc();
        }
        double pi =  calculator.getPI(peptide,true,true);
        return (Math.round(pi*100))/100.0;
    }
}

