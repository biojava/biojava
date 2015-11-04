package org.biojava.nbio.structure.io;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.BondImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompBond;
/**
 * 
 * Builds a list of connections from ligands in a Structure
 * and adds them to that Structure.
 * 
 * @author edlunde
 * @author larsonm
 * 
 * @since 10/30/2015.
 */
public class LigandConnectMaker {
	
	private Structure structure;
	
	public LigandConnectMaker( Structure struct ){
		structure = struct;
	}
	public void addLigandConnections() {

		ArrayList<Bond> bondsToAdd = new ArrayList<Bond>();
		
		for (Chain chain : structure.getChains()) {
			List<Group> groups = chain.getAtomLigands();

			for (Group group : groups) {
				// atoms with no residue number don't have atom information
				if (group.getResidueNumber() == null) {
					continue;
				}

				ChemComp compoundChemComp = ChemCompGroupFactory.getChemComp(group
						.getPDBName());
				
				for (ChemCompBond chemCompBond : compoundChemComp.getBonds()) {

					Atom a = group.getAtom(chemCompBond.getAtom_id_1());
					Atom b = group.getAtom(chemCompBond.getAtom_id_2());

					if ( a != null && b != null){
						int bondOrder = chemCompBond.getNumericalBondOrder();
						bondsToAdd.add(new BondImpl(a, b, bondOrder));
					} else  {
						// Some of the atoms were missing. That's fine, there's
						// nothing to do in this case.
					}
				}
			}
		}
		List<Map<String, Integer>> newBondMap = buildConectMap(bondsToAdd);
		List<Map<String, Integer>> existingConnections = structure.getConnections();
		existingConnections.addAll(newBondMap);
	}
	
	
	/**
	 * Starting with a list of all the bonds present in a model, create a CONECT record list
	 * that only includes the ligand connections of atoms present in the final structure.
	 * TODO : optionally include known hbonds.
	 * 
	 * @param bondList 
	 */
	private List<Map<String, Integer>> buildConectMap(List<Bond> bondList) {
		// Create a data structure of a HashMap to identify all the atom1 -> (set) connections first,
		// IBond list is a pairwise list only.
		HashMap<Integer, List<Integer>> bondMap = new HashMap<Integer, List<Integer>>();
		
		// Create a HashMap to show current atom indices present.
		HashMap<Integer, Group> atomMap = mapAtomIndices(structure);
		
		
		for (Bond bond : bondList) {
			final Atom atom1 = bond.getAtomA();
			final Atom atom2 = bond.getAtomB();
			
			final Group ag1 = atom1.getGroup();
			final Group ag2 = atom2.getGroup();
			
			// Save a bond to our map only if both atoms are ligands from the same group
			if ((atomMap.get(atom1.getPDBserial()) != null) && (atomMap.get(atom2.getPDBserial()) != null)
					&& ag1.equals(ag2) ){
			
				if (bondMap.get(atom1.getPDBserial()) == null) {
					bondMap.put(atom1.getPDBserial(), new ArrayList<Integer>());
				}
				// Add atom2 to the atom1 list.
				bondMap.get(atom1.getPDBserial()).add(atom2.getPDBserial());
			}
		}
		
		return buildBioJavaConnectFromMap(bondMap);
	}
	
	/**
	 * Convert a HashMap of ID: bonded ID1, ID2,.. IDN into BioJava CONECT map format.
	 */
	private List<Map<String, Integer>> buildBioJavaConnectFromMap(HashMap<Integer, List<Integer>> bondMap) {
		// restructure into the mapping used by BioJava
		List<Map<String, Integer>> connectList = new ArrayList<Map<String, Integer>>();
		/**
		 * the HashMap for a single CONECT line contains the following fields:
		 * atomserial (mandatory) : Atom serial number
		 * bond1 .. bond4 (optional): Serial number of bonded atom
		 * hydrogen1 .. hydrogen4 (optional):Serial number of hydrogen bonded atom
		 * salt1 .. salt2 (optional): Serial number of salt bridged atom
		 */
		
		// Create a HashMap per bond.  
		//  Add "atomserial" = atom1#
		//  Add "bond1" = atom2#
		//  Add "bond2", "bond3", "bond4" optionally
		//  Add "hydrogen1" - 4 optionally.
		//  Add "salt1" - 4 optionally.
		for (Integer atomA_index : bondMap.keySet()) {
			
			// Can create up to 4 bonds with same atomserial.
			List<Integer> bondedAtomList = bondMap.get(atomA_index);
			int freeBonds = 0;
			HashMap<String, Integer> connectRecord = null;
			boolean dirty = false;
			// Walk through the atom list.
			for (Integer atomB_index : bondedAtomList) {
				// When run out of space in a list, allocate a new list.
				if (freeBonds == 0) {
					// if we have a list.
					if (connectRecord != null) {
						connectList.add(connectRecord);
					}
					
					// Need another connect record.
					connectRecord = new HashMap<String, Integer>();
					dirty = false;
					freeBonds = 4;
					connectRecord.put("atomserial", atomA_index);
				}
				switch (freeBonds) {
					case 1: connectRecord.put("bond4", atomB_index); break;
					case 2: connectRecord.put("bond3", atomB_index); break;
					case 3: connectRecord.put("bond2", atomB_index); break;
					case 4: connectRecord.put("bond1", atomB_index); break;
				}
				
				dirty = true;
				freeBonds--;
			}
			// Need to add the last connectRecord
			if (dirty) {
				connectList.add(connectRecord);
			}
		}
		
		return connectList;
	}
	
	/**
	 * Build a lookup table of atom indices.
	 * @param index
	 * @return
	 */
	private HashMap<Integer, Group> mapAtomIndices(Structure structure) {
		HashMap<Integer, Group> atom_map = new HashMap<Integer, Group>();
		
		List<Chain> chains = structure.getChains();
		
		for (Chain a_chain : chains) {
			List<Group> groups = a_chain.getAtomGroups();
			for (Group a_group : groups) {
				List<Atom> atoms = a_group.getAtoms();
				for (Atom a_atom : atoms) {
					atom_map.put(a_atom.getPDBserial(), a_group);
				}
			}
		}
		
		return atom_map;
	}
}
