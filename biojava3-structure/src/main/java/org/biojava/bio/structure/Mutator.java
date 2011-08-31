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
 */

package org.biojava.bio.structure;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.structure.io.PDBParseException;


/** A class that can change one amino acid to another. Side chain atoms are neglected, only the Cb atom is kept.
 * 
 *
 * example usage:
 * <pre>
 String filename   =  "/Users/ap3/WORK/PDB/5pti.pdb" ;
 String outputfile =  "/Users/ap3/WORK/PDB/mutated.pdb" ;
 
 PDBFileReader pdbreader = new PDBFileReader();
 
 try{
     Structure struc = pdbreader.getStructure(filename);
     System.out.println(struc);
 
 
     String chainId = " ";
     String pdbResnum = "3";
     String newType = "ARG";
 
     // mutate the original structure and create a new one.
      Mutator m = new Mutator();
      Structure newstruc = m.mutate(struc,chainId,pdbResnum,newType);
  
      FileOutputStream out= new FileOutputStream(outputfile); 
      PrintStream p =  new PrintStream( out );
  
      p.println (newstruc.toPDB());
  
      p.close();
  
  
  } catch (Exception e) {
      e.printStackTrace();
  } 
  </pre>
  * @author Andreas Prlic
  * @since 1.5
  * @version %I% %G%
  */       
public class Mutator{
    List<String> supportedAtoms;
    
    public Mutator(){
        supportedAtoms = new ArrayList<String>();
        supportedAtoms.add("N");
        supportedAtoms.add("CA");
        supportedAtoms.add("C");
        supportedAtoms.add("O");
        supportedAtoms.add("CB");
    }
    
    /** creates a new structure which is identical with the original one. 
     * only one amino acid will be different.
     * 
     * 
     * 
     * 
     * @param struc the structure object that is the container for the residue to be mutated
     * @param chainId the id (name) of the chain to be mutated. @see Chain.getName()
     * @param pdbResnum the PDB residue number of the residue
     * @param newType the new residue type (3 characters)
     * @return a structure object where one residue has been modified
     * @throws PDBParseException
     */
    public Structure  mutate(Structure struc, String chainId, String pdbResnum, String newType) 
    throws PDBParseException{
        
        
        // create a  container for the new structure
        Structure newstruc = new StructureImpl();
        
        // first we need to find our corresponding chain
        
        // get the chains for model nr. 0
        // if structure is xray there will be only one "model".
        List<Chain> chains = struc.getChains(0);
        
        // iterate over all chains.
        Iterator<Chain> iter = chains.iterator();
        while (iter.hasNext()){
            Chain c = (Chain)iter.next();
            if (c.getChainID().equals(chainId)) {
                // here is our chain!
                
                Chain newchain = new ChainImpl();
                newchain.setChainID(c.getChainID());
                
                List<Group> groups = c.getAtomGroups();
                
                // now iterate over all groups in this chain.
                // in order to find the amino acid that has this pdbRenum.               
                
                Iterator<Group> giter = groups.iterator();
                while (giter.hasNext()){
                    Group g = (Group) giter.next();
                    String rnum = g.getResidueNumber().toString();
                    
                    // we only mutate amino acids
                    // and ignore hetatoms and nucleotides in this case                   
                    if (rnum.equals(pdbResnum) && (g.getType().equals("amino"))){
                        
                        // create the mutated amino acid and add it to our new chain
                        AminoAcid newgroup = mutateResidue((AminoAcid)g,newType);
                        newchain.addGroup(newgroup);
                    }
                    else {
                        // add the group  to the new chain unmodified.
                        newchain.addGroup(g);
                    }
                }
                
                // add the newly constructed chain to the structure;
                newstruc.addChain(newchain);
            } else {
                // this chain is not requested, add it to the new structure unmodified.
                newstruc.addChain(c);
            }
            
        }
        return newstruc;
    }
    
    /** create a new residue which is of the new type. 
     * Only the atoms N, Ca, C, O, Cb will be considered.
     * @param oldAmino
     * @param newType
     * @return a new, mutated, residue 
     * @throws PDBParseException
     */
    public AminoAcid mutateResidue(AminoAcid oldAmino, String newType)
    throws PDBParseException {
        
        AminoAcid newgroup = new AminoAcidImpl();
        
        newgroup.setResidueNumber(oldAmino.getResidueNumber());
        newgroup.setPDBName(newType);
        
        
        AtomIterator aiter =new AtomIterator(oldAmino);
        while (aiter.hasNext()){
            Atom a = (Atom)aiter.next();
            if ( supportedAtoms.contains(a.getName())){
                newgroup.addAtom(a);
            }
        }
        
        return newgroup;
        
    }
    
}
