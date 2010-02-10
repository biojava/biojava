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

import java.io.ObjectStreamException;
import java.io.Serializable;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;





/** The protease class stores parameters needed by Digest to digest a protein sequence.

 * A custom protease can be created or one derived from the attributes set in the

 * ProteaseManager.xml resource.

 * @author Michael Jones
 * @author Mark Schreiber (refactoring to ProteaseManager)

 */

public class Protease implements Serializable {

    public static final String TRYPSIN = ProteaseManager.TRYPSIN;
    public static final String LYS_C = ProteaseManager.LYS_C;
    public static final String ARG_C = ProteaseManager.ARG_C;
    public static final String ASP_N = ProteaseManager.ASP_N;
    public static final String GLU_C_BICARB = ProteaseManager.GLU_C_BICARB;
    public static final String GLU_C_PHOS = ProteaseManager.GLU_C_PHOS;
    public static final String CHYMOTRYP = ProteaseManager.CHYMOTRYP;
    public static final String CNBr = ProteaseManager.CNBr;

    private SymbolList cleavageResidues;

    private SymbolList notCleaveResidues;

    private boolean endoProtease = true;
    private String name;


    protected Protease(SymbolList cleaveRes,
                    boolean endoProtease,
                    SymbolList notCleaveRes,
                    String name)throws IllegalSymbolException{

      if(cleaveRes.getAlphabet() != ProteinTools.getAlphabet()){
        throw new IllegalSymbolException
        ("Cleaveage residues must be from the PROTEIN alphabet");
      }
      if(notCleaveRes.getAlphabet() != ProteinTools.getAlphabet() &&
         notCleaveRes.getAlphabet() != Alphabet.EMPTY_ALPHABET){
        throw new IllegalSymbolException
        ("Cleaveage residues must be from the PROTEIN alphabet or an EmptySymbolList");
      }
      this.cleavageResidues = cleaveRes;
      this.endoProtease = endoProtease;
      this.notCleaveResidues = notCleaveRes;
      this.name = name;
    }

    /**
     * @deprecated Creating a Protease with this constructor will not register it
     * with the ProteaseManager (use ProteaseManager.createProtease())
     */

    public Protease(SymbolList cleaveRes,

                    boolean endoProtease,

                    SymbolList notCleaveRes)

                                   throws IllegalSymbolException, BioException {



    }




    /**
     * @deprecated Creating a Protease with this constructor will not register it
     * with the ProteaseManager (use ProteaseManager.createProtease())
     */

    public Protease(String cleaveRes,

                    boolean endoProtease,

                    String notCleaveRes)

                                   throws IllegalSymbolException, BioException {

        this.cleavageResidues = createSymbolList(cleaveRes);

        this.endoProtease = endoProtease;

        this.notCleaveResidues = createSymbolList(notCleaveRes);

    }



    /**
     * @deprecated Creating a Protease with this constructor will not register it
     * with the ProteaseManager (use ProteaseManager.createProtease())
     */

    public Protease(String cleavageRes, boolean endoProtease)

                                  throws IllegalSymbolException, BioException {

        this.cleavageResidues = createSymbolList(cleavageRes);

        this.endoProtease = endoProtease;

        this.notCleaveResidues = createSymbolList("");

    }


    /**
     * The list of residues that the protease will cleave at.
     * @return the residues as a SymbolList
     */
    public SymbolList getCleaveageResidues()

    {
        return cleavageResidues;
    }

    /**
     * Gets the name of this Protease
     * @return the name as a String
     */
    public String getName(){
      return name;
    }


    /**
     * The list of residues that will prevent cleavage if they follow the cleavage
     * residue.
     */
    public SymbolList getNotCleaveResidues()

    {
        return notCleaveResidues;
    }



    public boolean isEndoProtease()

    {

        return endoProtease;

    }





    /**
     * Get the list of Protease names defined in the ProteaseManager
     * (Internally calls ProteaseManager.
     * @return A String array of protease names
     */

    public static String[] getProteaseList(){
      java.util.Set s = ProteaseManager.getNames();
      String[] names = new String[s.size()];
      return (String[])s.toArray(names);
    }



    /**
     * Retrieves a reference to the named Protease.
     * (Internally calls ProteaseManager.getProteaseByName())
     * @param proteaseName A protease name that is registered in the ProteaseManager (case sensitive)
     * @return A Protease instance for the given protease name
     */

    public static final Protease getProteaseByName(String proteaseName)
                                 throws BioException {
      return ProteaseManager.getProteaseByName(proteaseName);

    }



    private SymbolList createSymbolList(String seq)

                                  throws IllegalSymbolException, BioException {

        if(seq == null || seq.trim().equals("")){
          return SymbolList.EMPTY_LIST;
        }

        SymbolList sList;

        FiniteAlphabet prot

                 = (FiniteAlphabet)AlphabetManager.alphabetForName("PROTEIN");



        SymbolTokenization tokenization = prot.getTokenization("token");

        sList = new SimpleSymbolList (tokenization, seq);

        return sList;

    }

    /**
     * Prevent duplication of the object during Serialization
     * @throws ObjectStreamException
     */
    protected Object readResolve() throws ObjectStreamException{
      try {
        if(ProteaseManager.registered(this.getName())){
          return ProteaseManager.getProteaseByName(this.getName());
        }else{
          ProteaseManager.registerProtease(this);
          return this;
        }
      }
      catch (BioException ex) {
        throw new BioError(
          "Assertion Error: Cannot register Protease instance following de serialization", ex
        );
      }
    }

}

