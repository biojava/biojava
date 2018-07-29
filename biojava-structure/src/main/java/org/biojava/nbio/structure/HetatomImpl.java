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
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.GroupToSDF;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

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
public class HetatomImpl implements Group {

	private static final Logger logger = LoggerFactory.getLogger(HetatomImpl.class);

	private static final long serialVersionUID = 4491470432023820382L;

	/**
	 * The GroupType is HETATM
	 */
	public static final GroupType type = GroupType.HETATM ;

	private Map<String, Object> properties ;

	private long id;

	/** stores if 3d coordinates are available. */
	protected boolean pdb_flag ;

	/** 3 letter name of amino acid in pdb file. */
	protected String pdb_name ;

	protected ResidueNumber residueNumber;

	protected List<Atom> atoms ;

	private Chain parent;
	
	private boolean isHetAtomInFile;

	/**
	 * Behaviors for how to balance memory vs. performance.
	 * @author Andreas Prlic
	 */
	public static enum PerformanceBehavior {

		/** use a built-in HashMap for faster access to memory, at the price of more memory consumption */
		BETTER_PERFORMANCE_MORE_MEMORY,

		/** Try to minimize memory consumption, at the price of slower speed when accessing atoms by name */
		LESS_MEMORY_SLOWER_PERFORMANCE

	}

	public static PerformanceBehavior performanceBehavior=PerformanceBehavior.LESS_MEMORY_SLOWER_PERFORMANCE;

	private Map<String,Atom> atomNameLookup;

	protected ChemComp chemComp ;

	private List<Group> altLocs;

	/**
	 *  Construct a Hetatom instance.
	 */
	public HetatomImpl() {
		super();

		pdb_flag = false;
		pdb_name = null ;

		residueNumber = null;
		atoms    = new ArrayList<Atom>();
		properties = new HashMap<String,Object>();
		parent = null;
		chemComp = null;
		altLocs = null;

		if ( performanceBehavior == PerformanceBehavior.BETTER_PERFORMANCE_MORE_MEMORY)
			atomNameLookup = new HashMap<String,Atom>();
		else
			atomNameLookup = null;
	}


	/**
	 *  returns true or false, depending if this group has 3D coordinates or not.
	 * @return true if Group has 3D coordinates
	 */
	@Override
	public boolean has3D() {
		return pdb_flag;
	}

	/** flag if group has 3D data.
	 *
	 * @param flag  true to set flag that this Group has 3D coordinates
	 */
	@Override
	public void setPDBFlag(boolean flag){
		pdb_flag = flag ;
	}

	/** Set three character name of Group .
	 *
	 * @param s  a String specifying the PDBName value
	 * @see #getPDBName
	 */
	@Override
	public void setPDBName(String s) {
		// hetatoms can have pdb_name length < 3. e.g. CU (see 1a4a position 1200 )
		//if (s.length() != 3) {
		//throw new PDBParseException("amino acid name is not of length 3!");
		//}
		if (s != null && s.equals("?")) logger.info("invalid pdbname: ?");
		pdb_name =s ;

	}

	/**
	 * Returns the PDBName.
	 *
	 * @return a String representing the PDBName value
	 * @see #setPDBName
	 */
	@Override
	public String getPDBName() { return pdb_name;}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void addAtom(Atom atom){
		atom.setGroup(this);
		atoms.add(atom);
		// TODO this check is useless, coords are always !=null since they are initialized to 0,0,0 in AtomImpl constructor. We need to review this - JD 2016-09-14
		if (atom.getCoordsAsPoint3d() != null){
			// we have got coordinates!
			setPDBFlag(true);
		}

		if (atomNameLookup != null){

			Atom existingAtom = atomNameLookup.put(atom.getName(), atom);

			// if an atom with same name is added to the group that has to be some kind of problem,
			// we need to warn properly
			if (existingAtom != null) {
				String altLocStr = "";
				char altLoc = atom.getAltLoc();
				if (altLoc != ' ') altLocStr = "(alt loc '" + altLoc + "')";
				logger.warn("An atom with name " + atom.getName() + " " + altLocStr + " is already present in group: " + this.toString() + ". The atom with serial " + atom.getPDBserial() + " will be ignored in look-ups.");
			}
		}
	};


	/** remove all atoms
	 *
	 */
	@Override
	public void clearAtoms() {
		atoms.clear();
		setPDBFlag(false);
		if ( atomNameLookup != null)
			atomNameLookup.clear();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public int size(){ return atoms.size();   }

	/**
	 * {@inheritDoc}
	 */
	@Override
	public List<Atom> getAtoms(){
		return atoms ;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setAtoms(List<Atom> atoms) {

		// important we are resetting atoms to a new list, we need to reset the lookup too!
		if ( atomNameLookup != null)
			atomNameLookup.clear();

		for (Atom a: atoms){
			a.setGroup(this);
			if ( atomNameLookup != null)
				atomNameLookup.put(a.getName(),a);
		}
		this.atoms = atoms;
		if (!atoms.isEmpty()) {
			pdb_flag = true;
		}

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Atom getAtom(String name) {
		if ( atomNameLookup != null)
			return atomNameLookup.get(name);
		else {
			/** This is the performance penalty we pay for NOT using the atomnameLookup in PerformanceBehaviour.LESS_MEMORY_SLOWER_PERFORMANCE
			 */
			for (Atom a : atoms) {
				if (a.getName().equals(name)) {
					return a;
				}
			}
			return null;
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Atom getAtom(int position) {

		if ((position < 0)|| ( position >= atoms.size())) {
			//throw new StructureException("No atom found at position "+position);
			return null;
		}
		return atoms.get(position);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean hasAtom(String fullName) {

		if ( atomNameLookup != null) {
			Atom a = atomNameLookup.get(fullName.trim());
			return a != null;
		} else {
			/** This is the performance penalty we pay for NOT using the atomnameLookup in PerformanceBehaviour.LESS_MEMORY_SLOWER_PERFORMANCE
			 */
			for (Atom a : atoms) {
				if (a.getName().equals(fullName)) {
					return true;
				}
			}
			return false;


		}

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public GroupType getType(){ return type;}

	@Override
	public String toString(){

		String str = "Hetatom "+ residueNumber + " " + pdb_name +  " "+ pdb_flag;
		if (pdb_flag) {
			str = str + " atoms: "+atoms.size();
		}
		if ( altLocs != null)
			str += " has altLocs :" + altLocs.size();


		return str ;

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean hasAminoAtoms(){
		// if this method call is performed too often, it should become a
		// private method and provide a flag for Group object ...

		return hasAtom(StructureTools.CA_ATOM_NAME) &&
				hasAtom(StructureTools.C_ATOM_NAME) &&
				hasAtom(StructureTools.N_ATOM_NAME) &&
				hasAtom(StructureTools.O_ATOM_NAME);

	}

	@Override
	public boolean isPolymeric() {

		ChemComp cc = getChemComp();

		if ( cc == null)
			return getType().equals(GroupType.AMINOACID) || getType().equals(GroupType.NUCLEOTIDE);

		ResidueType rt = cc.getResidueType();

		if ( rt.equals(ResidueType.nonPolymer))
			return false;

		PolymerType pt = rt.getPolymerType();

		return PolymerType.PROTEIN_ONLY.contains(pt) ||
				PolymerType.POLYNUCLEOTIDE_ONLY.contains(pt) ||
				ResidueType.lPeptideLinking.equals(rt);


	}

	@Override
	public boolean isAminoAcid() {

		ChemComp cc = getChemComp();

		if ( cc == null)
			return getType().equals(GroupType.AMINOACID);


		ResidueType rt = cc.getResidueType();

		if ( rt.equals(ResidueType.nonPolymer))
			return false;

		PolymerType pt = rt.getPolymerType();

		return PolymerType.PROTEIN_ONLY.contains(pt);

	}

	@Override
	public boolean isNucleotide() {

		ChemComp cc = getChemComp();

		if ( cc == null)
			return  getType().equals(GroupType.NUCLEOTIDE);

		ResidueType rt = cc.getResidueType();

		if ( rt.equals(ResidueType.nonPolymer))
			return false;

		PolymerType pt = rt.getPolymerType();

		return PolymerType.POLYNUCLEOTIDE_ONLY.contains(pt);


	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setProperties(Map<String,Object> props) {
		properties =  props ;
	}

	/** return properties.
	 *
	 * @return a HashMap object representing the properties value
	 * @see #setProperties
	 */
	@Override
	public Map<String, Object> getProperties() {
		return properties ;
	}

	/** set a single property .
	 *
	 * @see #getProperties
	 * @see #getProperty
	 */
	@Override
	public void setProperty(String key, Object value){
		properties.put(key,value);
	}

	/** get a single property .
	 * @param key  a String
	 * @return an Object
	 * @see #setProperty
	 * @see #setProperties
	 */
	@Override
	public Object getProperty(String key){
		return properties.get(key);
	}


	/** return an AtomIterator.
	 *
	 * @return an Iterator object
	 */
	@Override
	public Iterator<Atom> iterator() {
		return new AtomIterator(this);
	}

	/** returns and identical copy of this Group object .
	 * @return  and identical copy of this Group object
	 */
	@Override
	public Object clone() {

		HetatomImpl n = new HetatomImpl();
		n.setPDBFlag(has3D());
		n.setResidueNumber(residueNumber);

		n.setPDBName(getPDBName());

		//clone atoms and bonds.
		cloneAtomsAndBonds(n);
		
		// copying the alt loc groups if present, otherwise they stay null
		if (altLocs!=null) {
			for (Group altLocGroup:this.altLocs) {
				Group nAltLocGroup = (Group)altLocGroup.clone();
				n.addAltLoc(nAltLocGroup);
			}
		}
		
		if (chemComp!=null)
			n.setChemComp(chemComp);

		return n;
	}


	protected void cloneAtomsAndBonds(Group newGroup) {
		// copy the atoms
		for (Atom atom1 : atoms) {
			Atom atom = (Atom) atom1.clone();
			newGroup.addAtom(atom);
			atom.setGroup(newGroup);
		}
		// copy the bonds
		for (int i=0;i<atoms.size();i++) {
			Atom atom1 = atoms.get(i);
			List<Bond> bonds1 = atom1.getBonds();
			if (bonds1 != null) {
				for (Bond b : bonds1) {
					int atomAIndex = atoms.indexOf(b.getAtomA());
					int atomBIndex = atoms.indexOf(b.getAtomB());
					// The order of the atoms are the same on the original and the cloned object, which we use here.
					Bond newBond = new BondImpl(newGroup.getAtom(atomAIndex), newGroup.getAtom(atomBIndex), b.getBondOrder(), false);
					newGroup.getAtom(i).addBond(newBond);
				}
			}
		}
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

	@Override
	public ChemComp getChemComp() {
		if  ( chemComp == null ) {
			chemComp = ChemCompGroupFactory.getChemComp(pdb_name);
			if (chemComp == null) logger.info("getChemComp: " + pdb_name);
		}
		return chemComp;
	}

	@Override
	public void setChemComp(ChemComp cc) {
		chemComp = cc;

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public void setChain(Chain chain) {
		this.parent = chain;
		//TODO: setChain(), getChainName() and ResidueNumber.set/getChainName() are
		//duplicating functionality at present and could give different values.
		if (residueNumber != null) {
			residueNumber.setChainName(chain.getName());
		}

	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public Chain getChain() {
		return parent;
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public String getChainId() {
		if (parent == null) {
			return "";
		}
		return parent.getId();
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public ResidueNumber getResidueNumber() {

		return residueNumber;
	}


	@Override
	public void setResidueNumber(ResidueNumber residueNumber) {
		this.residueNumber = residueNumber;
	}

	@Override
	public void setResidueNumber(String chainId, Integer resNum, Character iCode) {
		this.residueNumber = new ResidueNumber(chainId, resNum, iCode);
	}

	@Override
	public boolean hasAltLoc() {
		if ( altLocs == null)
			return false;
		return !altLocs.isEmpty();
	}

	@Override
	public List<Group> getAltLocs() {
		if ( altLocs == null)
			return new ArrayList<Group>();
		return altLocs;
	}

	@Override
	public Group getAltLocGroup(Character altLoc) {

		Atom a = getAtom(0);
		if ( a == null) {
			return null;
		}

		// maybe the alt loc group in question is myself
		if (a.getAltLoc().equals(altLoc)) {
			return this;
		}

		if (altLocs == null || altLocs.isEmpty())
			return null;

		for (Group group : altLocs) {
			if (group.getAtoms().isEmpty())
				continue;

			// determine this group's alt-loc character code by looking
			// at its first atom's alt-loc character
			Atom b = group.getAtom(0);
			if ( b == null)
				continue;

			if (b.getAltLoc().equals(altLoc)) {
				return group;
			}
		}

		return null;
	}

	@Override
	public void addAltLoc(Group group) {
		if ( altLocs == null) {
			altLocs = new ArrayList<Group>();
		}
		altLocs.add(group);

	}

	@Override
	public boolean isWater() {
		return GroupType.WATERNAMES.contains(pdb_name);
	}

	/** attempts to reduce the memory imprint of this group by trimming
	 * all internal Collection objects to the required size.
	 *
	 */
	@Override
	public void trimToSize(){

		if ( atoms instanceof ArrayList<?>) {
			ArrayList<Atom> myatoms = (ArrayList<Atom>) atoms;
			myatoms.trimToSize();
		}
		if ( altLocs instanceof ArrayList<?>){
			ArrayList<Group> myAltLocs = (ArrayList<Group>) altLocs;
			myAltLocs.trimToSize();
		}

		if ( hasAltLoc()) {
			for (Group alt : getAltLocs()){
				alt.trimToSize();
			}
		}

		// now let's fit the hashmaps to size
		properties = new HashMap<String, Object>(properties);

		if ( atomNameLookup != null)
			atomNameLookup = new HashMap<String,Atom>(atomNameLookup);

	}


	@Override
	public String toSDF() {
		// Function to return the SDF of a given strucutre
		GroupToSDF gts = new GroupToSDF();
		return gts.getText(this);
	}

	@Override
	public boolean isHetAtomInFile() {
		return isHetAtomInFile;
	}
	
	@Override
	public void setHetAtomInFile(boolean isHetAtomInFile) {
		this.isHetAtomInFile = isHetAtomInFile;
	}




}
