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
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolPropertyTable;

/**
 * <code>MassCalc</code> calculates the mass of peptides which for our
 * purposes are <code>SymbolList</code>s which contain
 * <code>Symbol</code>sfrom the protein <code>Alphabet</code>. It uses
 * the mono-isotopic and average-isotopic masses identical to those
 * specified at www.micromass.co.uk
 * Note: This class does not handle selenocysteine and pyrrolysine.
 *
 * @author M. Jones sdfsd
 * @author Keith James (minor changes)
 * @author Mark Schreiber
 * @author George Waldon - getMolecularWeight
 */
public class MassCalc {
    //Possibly these should be put in a configurable table somewhere
    /**
     * Constant value of Carbon monoisotopic mass
     */
    public static final double Cmono = 12.00000;

    /**
     * Constant value of  Hydrogen monoisotopic mass
     */
    public static final double Hmono = 1.0078250;

    /**
     * Constant value of Nitrogen monoisotopic mass
     */
    public static final double Nmono = 14.0030740;

    /**
     * Constant value of Oxygen monoisotopic mass
     */
    public static final double Omono = 15.9949146;

    /**
     * Constant value of Carbon average mass
     */
    public static final double Cavg = 12.011;

    /**
     * Constant value of Hydrogen average mass
     */
    public static final double Havg = 1.00794;

    /**
     * Constant value of Nitrogen average mass
     */
    public static final double Navg = 14.00674;

    /**
     * Constant value of Oxygen average mass
     */
    public static final double Oavg = 15.9994;

    //Save values here so that modifications are not global
    private HashMap mSymbolPropertyHash;


    private HashMap mVariableModPropertyHash;

    /*
     * If instance methods are being used this will store the
     * isotopically and MH_PLUS correct terminal mass to be created.
     * Need to be able to handle N and C term mods in the future.
     */
     private double termMass;

    /*
     * Creates new MassCalc.
     * @param isotopicType The type of isotopes to calculate. Either
     * mono isotopic or average isotopic. It realy is the name of
     * SymbolProperty table. Ex. SymbolPropertyTable.AVG_MASS or
     * SymbolPropertyTable.MONO_MASS */


    /**
     * Creates a new <code>MassCalc</code>.
     *
     * @param isotopicType a <code>String</code>. The type of isotopes
     * to calculate. Either mono isotopic or average
     * isotopic. Acceptable values are
     * <code>SymbolPropertyTable.AVG_MASS</code> or
     * <code>SymbolPropertyTable.MONO_MASS</code>.
     * @param MH_PLUS a <code>boolean</code>.
     */
    public MassCalc(String isotopicType, boolean MH_PLUS) {
        //Calculate hydroxyl mass
        termMass = calcTermMass(isotopicType, MH_PLUS);
        mSymbolPropertyHash = new HashMap();

        SymbolPropertyTable symbolPropertyTable =
            ProteinTools.getSymbolPropertyTable(isotopicType);

        mVariableModPropertyHash = new HashMap();

        Iterator symbolList = ProteinTools.getAlphabet().iterator();

        for(; symbolList.hasNext(); )
        {
            Symbol sym = (Symbol) symbolList.next();
            //SELENOCYSTEINE and PYRROLYSINE
            if(sym==ProteinTools.o() || sym==ProteinTools.u())
                continue;
            try {
                try {
                    Double value =
                        new Double(symbolPropertyTable.getDoubleValue(sym));
                    mSymbolPropertyHash.put(sym, value);
                } catch (NullPointerException npe) {
                    //This seems to be happening when a amino acid is
                    //in the property table that doesn't have a residue value
                }
            }
            catch (IllegalSymbolException ise)
            {
                throw new BioError("Error setting properties for Symbol " + sym, ise);
            }
        }
    }

    /**
     * Use this to set a post translational modification for the
     * <code>Symbol</code> represented by this character. It will only
     * affect the current <code>MassCalc</code> instance and will not
     * affect the static method.
     *
     * @param symbolToken a <code>char</code> representing a
     * <code>Symbol</code>.
     * @param mass a <code>double</code> to be the new mass of the
     * residue.
     *
     * @exception IllegalSymbolException if the <code>Symbol</code> is
     * not recognised.
     */
    public void setSymbolModification(char symbolToken, double mass)
        throws IllegalSymbolException
    {
        SymbolTokenization toke;

        try {
            toke = ProteinTools.getAlphabet().getTokenization("token");
        } catch (BioException ex) {
            throw new BioError("Expected a tokenization", ex);
        }

        Symbol sym = toke.parseToken("" + symbolToken);
        mSymbolPropertyHash.put(sym, new Double(mass));
    }

    /** Add Variable modifications.  If multiple masses are
     * set by this method more then one mass will be returned for a mass
     * calculation. For example if a peptide contains two Mets and the user sets
     * the native and oxidized mass for the Met then the masses returned will be
     * of the peptide with 0, 1 and 2 modified Mets.
     * @param residue The one char id for this residue
     * @param masses
     * @throws IllegalSymbolException
     * @see #getVariableMasses(SymbolList peptide)
     * @see #addVariableModification(Symbol residue,double[] masses)
     */
     public void addVariableModification(char residue,
                                         double[] masses)
        throws IllegalSymbolException{
        Symbol sym = getSymbolForChar(residue);
        addVariableModification(sym, masses);
     }


     /** Add Variable modifications. If multiple masses are
      * set by this method more then one mass will be returned for a mass
      * calculation. For example if a peptide contains two Mets and the user sets
      * the native and oxidized mass for the Met then the masses returned will be
      * of the peptide with 0, 1 and 2 modified Mets.
      */
     public void addVariableModification(Symbol residue,
                                         double[] masses)
        throws IllegalSymbolException{

        List massList = new LinkedList();
        for(int i=0; i<masses.length; i++){
            massList.add(new Double(masses[i]));
        }
        getVariableModMap().put(residue, massList);
     }

     /**
      * Remove all variable modifications assocaited with this residue.
      *
      */
     public boolean removeVariableModifications(char residue)
        throws IllegalSymbolException{
        Symbol sym = getSymbolForChar(residue);
        return removeVariableModifications(sym);
     }

     /**
      * Remove all variable modifications assocaited with this residue.
      *
      */
     public boolean removeVariableModifications(Symbol residue){
        if(getVariableModMap().remove(residue) != null){
            return true;
        }else{
            return false;
        }
     }

     /** Calculate the molecular weight of a protein, making estimates whenever it is 
      * possible like averaging mass values for ambiguity symbols or counting
      * zero when gaps are encountered. 
      * The method is tolerant for ambiguity symbols as long as they can be
      * resolved to a series of atomic symbols whose mass is available in the 
      * ResidueProperties.xml configuration file or they are gaps.
      * The method returns the same value as getMass(SymbolList proteinSeq,
      * SymbolPropertyTable.AVG_MASS, false) when only atomic symbols are found
      * in the polypeptide.
      *
      * @since 1.5
      */
     public static final double getMolecularWeight(SymbolList proteinSeq)
     throws IllegalSymbolException {
        double pepMass = 0.0;
        Symbol gap = AlphabetManager.alphabetForName("PROTEIN").getGapSymbol();
        SymbolPropertyTable sPT =
            ProteinTools.getSymbolPropertyTable(SymbolPropertyTable.AVG_MASS);

        for (Iterator it = proteinSeq.iterator(); it.hasNext(); ) {
            Symbol s = (Symbol) it.next();
            if( s instanceof AtomicSymbol) {
                pepMass += sPT.getDoubleValue(s);
            } else {
                FiniteAlphabet matches = (FiniteAlphabet) s.getMatches();
                if(matches.size() == 0) {
                    if(s.equals(gap)) {
                        continue;
                    }
                    throw new IllegalSymbolException(
                        "Symbol " + s.getName() + " has no mass associated");
                } else {
                    int count = 0;
                    double mass = 0.0;
                    for(Iterator i = matches.iterator(); i.hasNext(); ) {
                        AtomicSymbol as = (AtomicSymbol) i.next();
                        //SELENOCYSTEINE and PYRROLYSINE
                        if(as==ProteinTools.o() ||
                                as==ProteinTools.u())
                            continue;
                        mass += sPT.getDoubleValue(as);
                        count++;
                    }
                    pepMass += mass/((double)count);
                }
            }
        }

        //Calculate hydroxyl mass
        double termMass = calcTermMass(SymbolPropertyTable.AVG_MASS, false);

        if (pepMass != 0.0) {
            pepMass += termMass;
        }

        return pepMass;
     }
     
    /**
     * <code>getMass</code> calculates the mass of this peptide. This
     * only works for the values in the ResidueProperties.xml
     * configuration file. It is probably slightly faster than the
     * instance method, but it does not handle post-translational
     * modifications.
     *
     * @param proteinSeq a <code>SymbolList</code> whose mass is to be
     * calculated. This should use the protein alphabet.
     * @param isotopicType a <code>String</code> The type of isotopes
     * to calculate. Either mono isotopic or average
     * isotopic. Acceptable values are
     * <code>SymbolPropertyTable.AVG_MASS</code> or
     * <code>SymbolPropertyTable.MONO_MASS</code>.
     * @param MH_PLUS a <code>boolean</code> true if the value needed
     * is the MH+ mass.
     *
     * @return a <code>double</code> mass of the peptide.
     *
     * @exception IllegalSymbolException if the
     * <code>SymbolList</code> contains illegal
     * <code>Symbol</code>s.
     */
     public static final double getMass(SymbolList proteinSeq,
             String isotopicType,
             boolean MH_PLUS)
        throws IllegalSymbolException
    {

        double pepMass = 0.0;

        SymbolPropertyTable sPT =
            ProteinTools.getSymbolPropertyTable(isotopicType);

        for (Iterator it = proteinSeq.iterator(); it.hasNext(); ) {
            Symbol s = (Symbol) it.next();
            if(!(s instanceof AtomicSymbol)) {
                throw new IllegalSymbolException(
                    "Symbol " + s.getName() + " is not atomic");
            }
            pepMass += sPT.getDoubleValue(s);
        }

        //Calculate hydroxyl mass
        double termMass = calcTermMass(isotopicType, MH_PLUS);

        if (pepMass != 0.0) {
            pepMass += termMass;
        }

        return pepMass;
    }

    /**
     * Get the Mass of this peptide. Use this if you want to set fixed
     * modifications and have created an instance of MassCalc. The
     * value is calculated using the value of MH_PLUS defined in the
     * constructor. The static method may be faster.
     *
     * @param proteinSeq The sequence for mass calculation
     *
     * @return The mass of the sequence */
    public double getMass(SymbolList proteinSeq)
        throws IllegalSymbolException
    {

        double pepMass = 0.0;

        HashMap symbolPropertyMap = getSymbolPropertyMap();

        for (Iterator it = proteinSeq.iterator(); it.hasNext(); ) {
            Symbol s = (Symbol) it.next();
            if(! symbolPropertyMap.containsKey(s)){
                throw new IllegalSymbolException(s, "The mass of the symbol "+s.getName()+" is unknown");
            }
            Double mass = (Double) symbolPropertyMap.get(s);
            pepMass += mass.doubleValue();
        }
        pepMass += getTermMass();

        return pepMass;
    }

     /**
      * Get all masses including the variable mass.
      * Allgorythm
      *
      * 1 Get the first residue of the sequence
      * create a list of all the standard and non-standard massses for this reidue
      * for each residue mass goto 1 with the sequence of all residues after the current residue
      *     add the residue mass to each mass from 1 to the list
      *
      *
      *     @see #addVariableModification
      */
     public double[] getVariableMasses(SymbolList peptide) throws IllegalSymbolException {
        double[] vMasses = getVMasses(peptide);
        for(int i=0; i<vMasses.length; i++){
            vMasses[i] +=  getTermMass();
        }

        return vMasses;
     }


     private HashMap getVariableModMap() {
         return mVariableModPropertyHash;
     }


    private HashMap getSymbolPropertyMap(){
        return mSymbolPropertyHash;// = ProteinTools.getSymbolPropertyTable(name);
    }

    /**
     * <code>getTermMass</code> returns the terminal mass being used
     * by the instance methods.
     *
     * @return a <code>double</code> mass.
     */
    public double getTermMass(){
        return termMass;
    }

  /**
      *
      */
     private double[] getVMasses(SymbolList peptide) throws IllegalSymbolException {
        Set allMassList = new HashSet();

        Symbol sym = peptide.symbolAt(1);
        if(!getSymbolPropertyMap().containsKey(sym)){
            String msg = "No mass Set for Symbol " + sym.getName();
            throw new IllegalSymbolException(msg);
        }

        //Create a list for all masses of the current residue
        List curResMasses = null;
        if(getVariableModMap().containsKey(sym)){
            curResMasses = new LinkedList((List) getVariableModMap().get(sym));
        }else{
            curResMasses = new LinkedList();
        }
        curResMasses.add(getSymbolPropertyMap().get(sym));

        //Move through all masses and calculate the masses of all of the sub peptides
        Iterator it = curResMasses.iterator();
        while(it.hasNext()){
            double resMass = ((Double)it.next()).doubleValue();
            if(peptide.length() == 1){
                allMassList.add(new Double(resMass));
            }else{
                //Get all masses of remaining peptide
                double[] subMasses = getVMasses(peptide.subList(2, peptide.length()));
                //Get next modified mass for symbol
                for(int i=0; i<subMasses.length; i++){
                    double pepMass = resMass + subMasses[i];
                    allMassList.add(new Double(pepMass));
                }
            }
        }
        //Convert list to an array
        double masses[] = new double[allMassList.size()];
        int i=0;
        for(Iterator mit = allMassList.iterator(); mit.hasNext(); i++){
            masses[i] = ((Double)mit.next()).doubleValue();
        }
        return masses;
     }


    private static double calcTermMass(String isotopicType, boolean MH_PLUS) {
        double termMass = 0.0;
        if (isotopicType.equals(SymbolPropertyTable.AVG_MASS)) {
            //Add the C-terminal OH and N-Term H
            termMass += Havg + Oavg + Havg;
            //Add the extra H
            if (MH_PLUS) {
                termMass += Havg;
            }
        }
        else if (isotopicType.equals(SymbolPropertyTable.MONO_MASS)) {
            //Add the C-terminal OH and N-Term H
            termMass += Hmono + Omono + Hmono;
            //Add the extra H
            if (MH_PLUS) {
                termMass += Hmono;
            }
        }
        return termMass;
    }

     private Symbol getSymbolForChar(char symbolToken)
         throws IllegalSymbolException{
         SymbolTokenization toke;
         try {
             toke = ProteinTools.getAlphabet().getTokenization("token");
         } catch (BioException ex) {
             throw new BioError("Expected a tokenization", ex);
         }

         Symbol sym = toke.parseToken("" + symbolToken);
        return sym;
     }
}
