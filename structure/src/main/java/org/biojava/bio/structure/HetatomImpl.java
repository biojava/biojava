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
 * Created on 05.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.io.PDBParseException;

/**
 *
 * Generic Implementation of a Group interface.
 * AminoAcidImpl and NucleotideImpl are closely related classes.
 * @see AminoAcidImpl
 * @see NucleotideImpl
 * @author Andreas Prlic
 * @author Horvath Tamas
 * @version %I% %G%
 * @since 1.4
 */
public class HetatomImpl implements Group,Serializable {

   /**
    *
    */
   private static final long serialVersionUID = 4491470432023820382L;

   /** this is a "hetatm".
    *
    */
   public static final String type = GroupType.HETATM ;

   Map<String, Object> properties ;

   long id;

   /* stores if 3d coordinates are available. */
   protected boolean pdb_flag ;

   /* 3 letter name of amino acid in pdb file. */
   protected String pdb_name ;

   /* pdb numbering. */
   protected String pdb_code ;

   protected List<Atom> atoms ;

   Chain parent;

   Map<String,Atom> atomLookup = new HashMap<String,Atom>();

   /* Construct a Hetatom instance. */
   public HetatomImpl() {
      super();


      pdb_flag = false;
      pdb_name = null ;
      pdb_code = null ;
      atoms    = new ArrayList<Atom>();
      properties = new HashMap<String,Object>();
      parent = null;
   }

   /* returns an identical copy of this structure
     public Object clone() {
     Hetatom n = new Hetatom();
     }
    */




   /**
    *  returns true or false, depending if this group has 3D coordinates or not.
    * @return true if Group has 3D coordinates
    */
   public boolean has3D() {
      return pdb_flag;
   }

   /** flag if group has 3D data.
    *
    * @param flag  true to set flag that this Group has 3D coordinates
    */
   public void setPDBFlag(boolean flag){
      pdb_flag = flag ;
   }

   /**
    * Returns the PDBCode.
    * @see #setPDBCode
    * @return a String representing the PDBCode value
    */
   public String getPDBCode() {
      return pdb_code;
   }

   /** set the PDB code.
    * @see #getPDBCode
    */
   public void setPDBCode(String pdb) {
      pdb_code = pdb ;
   }

   /** set three character name of Group .
    *
    * @param s  a String specifying the PDBName value
    * @see #getPDBName
    * @throws PDBParseException ...
    */
   public void setPDBName(String s)
   throws PDBParseException
   {
      // hetatoms can have pdb_name length < 3. e.g. CU (see 1a4a position 1200 )
      //if (s.length() != 3) {
      //throw new PDBParseException("amino acid name is not of length 3!");
      //}
      pdb_name =s ;
   }

   /**
    * Returns the PDBName.
    *
    * @return a String representing the PDBName value
    * @see #setPDBName
    */
   public String getPDBName() { return pdb_name;}

   /** add an atom to this group. */
   public void addAtom(Atom atom){
      atom.setParent(this);
      atoms.add(atom);
      if (atom.getCoords() != null){
         // we have got coordinates!
         setPDBFlag(true);
      }
      atomLookup.put(atom.getName(),atom);
   };


   /** remove all atoms
    *
    */
   public void clearAtoms() {
      atoms.clear();
      setPDBFlag(false);
      atomLookup.clear();
   }

   /** getnumber of atoms.
    *  @return number of atoms
    */
   public int size(){ return atoms.size();   }

   /** get all atoms of this group .
    * returns a List of all atoms in this Group
    * @return an List object representing the atoms value
    */
   public List<Atom> getAtoms(){
      //Atom[] atms = (Atom[])atoms.toArray(new Atom[atoms.size()]);

      return atoms ;
   }


   /** set the atoms of this group
    * @see org.biojava.bio.structure.Atom
    * @param atoms a list of atoms
    */
   public void setAtoms(List<Atom> atoms){
      for (Atom a: atoms){
         a.setParent(this);
         atomLookup.put(a.getName(), a);
      }
      this.atoms = atoms;
      if ( atoms.size() > 0) {
         pdb_flag = true;
      }

   }


   /**  get an atom throws StructureException if atom not found.
    * @param name  a String
    * @return an Atom object
    * @throws StructureException ...
    */
   public Atom getAtom(String name)
   throws StructureException
   {
      // todo: add speedup by internal hashmap...

      Atom a = atomLookup.get(name);
      if ( a != null)
         return a;

      // something is deeply wrong...


      for (int i=0;i<atoms.size();i++){
         Atom atom = atoms.get(i);

         if (atom.getName().equals(name)){
            if ( name.equals("CA")) {
               if (atom.getElement().equals(Element.C))
               return atom;
            }   
         }

      }

      throw new StructureException(" No atom "+name + " in group " + pdb_name + " " + pdb_code + " !");

   }

   /**  Get an atom by the full PDB name e.g. " N  " for N. Throws StructureException if atom not found.
    * @param name  a String
    * @return an Atom object
    * @throws StructureException ...
    */
   public Atom getAtomByPDBname(String name)
   throws StructureException
   {

      for (int i=0;i<atoms.size();i++){
         Atom atom = atoms.get(i);
         if (atom.getFullName().equals(name)){
            return atom;
         }
      }

      throw new StructureException(" No atom "+name + " in group " + pdb_name + " " + pdb_code + " !");

   }

   /** return an atom by its position in the internal List.
    *
    * @param position  an int
    * @return an Atom object
    * @throws StructureException ...
    */
   public Atom getAtom(int position)
   throws StructureException
   {
      if ((position < 0)|| ( position >= atoms.size())) {
         throw new StructureException("No atom found at position "+position);
      }
      Atom a = atoms.get(position);
      return a ;
   }

   /** test is an Atom with name is existing. */
   public boolean hasAtom(String name){

      Atom a = atomLookup.get(name);
      if ( a != null)
         return true;

      return false;

      //       for (int i=0;i<atoms.size();i++){
      //            Atom atom = atoms.get(i);
      //            if (atom.getName().equals(name)){
      //                return true;
      //            }
      //        }
      //        return false ;
   }

   /**
    * Returns the type value.
    *
    * @return a String representing the type value
    */
   public String getType(){ return type;}

   public String toString(){

      String str = "Hetatom "+ pdb_code + " " + pdb_name +  " "+ pdb_flag;
      if (pdb_flag) {
         str = str + " atoms: "+atoms.size();
      }


      return str ;

   }



   /** calculate if a groups has all atoms required for an amino acid
     this allows to include chemically modified amino acids that
     are labeled hetatoms into some computations ... the usual way
     to identify if a group is an amino acid is getType() !
     <p>
     amino atoms are : N, CA, C, O, CB
     GLY does not have CB (unless we would calculate some artificially
     </p>

     Example: 1DW9 chain A first group is a Selenomethionine, provided as HETATM, but here returns true.
     <pre>
     HETATM    1  N   MSE A   1      11.720  20.973   1.584  0.00  0.00           N
     HETATM    2  CA  MSE A   1      10.381  20.548   1.139  0.00  0.00           C
     HETATM    3  C   MSE A   1       9.637  20.037   2.398  0.00  0.00           C
     HETATM    4  O   MSE A   1      10.198  19.156   2.985  0.00  0.00           O
     HETATM    5  CB  MSE A   1      10.407  19.441   0.088  0.00  0.00           C
     </pre>
     @see #getType

    */


   public boolean hasAminoAtoms(){
      // if this method call is performed too often, it should become a
      // private method and provide a flag for Group object ...

      String[] atoms ;
      if ( getType().equals("amino") & getPDBName().equals("GLY")){
         atoms = new String[] { "N","CA","C","O"};
      } else {
         atoms = new String[] { "N","CA","C","O","CB" };
      }


      for (int i = 0 ; i < atoms.length; i++) {
         if ( ! hasAtom(atoms[i])) {
            //System.out.println("not amino atoms");
            return false ;
         }
      }

      return true ;
   }


   /** properties of this amino acid. currerntly available properties.
    * are:
    * phi
    * psi
    *
    * @see #getProperties
    */
   public void setProperties(Map<String,Object> props) {
      properties =  props ;
   }

   /** return properties.
    *
    * @return a HashMap object representing the properties value
    * @see #setProperties
    */
   public Map<String, Object> getProperties() {
      return properties ;
   }

   /** set a single property .
    *
    * @see #getProperties
    * @see #getProperty
    */
   public void setProperty(String key, Object value){
      properties.put(key,value);
   }

   /** get a single property .
    * @param key  a String
    * @return an Object
    * @see #setProperty
    * @see #setProperties
    */
   public Object getProperty(String key){
      return properties.get(key);
   }


   /** return an AtomIterator.
    *
    * @return an Iterator object
    */
   public Iterator<Atom> iterator() {
      Iterator<Atom> iter = new AtomIterator(this);
      return iter ;
   }

   /** returns and identical copy of this Group object .
    * @return  and identical copy of this Group object
    */
   public Object clone(){

      HetatomImpl n = new HetatomImpl();
      n.setPDBFlag(has3D());
      n.setPDBCode(getPDBCode());
      try {
         n.setPDBName(getPDBName());
      } catch (PDBParseException e) {
         e.printStackTrace();
      }
      // copy the atoms
      for (int i=0;i<atoms.size();i++){
         Atom atom = atoms.get(i);
         n.addAtom((Atom)atom.clone());
      }
      return n;
   }

   /** Set the back-reference (to its parent Chain)
    * @param parent the parent Chain
    */

   public void setParent(Chain parent) {
      this.parent = parent ;
   }

   /** Returns the parent Chain of the Group
    *
    * @return Chain the Chain object that contains the Group
    *
    *
    */

   public Chain getParent() {
      return parent;
   }

   /** the Hibernate database ID
    *
    * @return the id
    */
   public long getId() {
      return id;
   }

   /** the Hibernate database ID
    *
    * @param id the hibernate id
    */
   public void setId(long id) {
      this.id = id;
   }



}
