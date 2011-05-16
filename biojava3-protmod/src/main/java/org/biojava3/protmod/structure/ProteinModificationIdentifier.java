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
 * Created on Jun 12, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.structure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;


import org.biojava3.protmod.Component;
import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModificationLinkage;
import org.biojava3.protmod.ProteinModification;
import org.biojava3.protmod.ProteinModificationRegistry;

/**
 * Identify attachment modification in a 3-D structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ProteinModificationIdentifier {
	

	private double bondLengthTolerance ;
	private boolean recordUnidentifiableModifiedCompounds ;
	private boolean recordAdditionalAttachments ;
	
	private Set<ModifiedCompound> identifiedModifiedCompounds = null;
	private Set<StructureAtomLinkage> unidentifiableAtomLinkages = null;
	private Set<StructureGroup> unidentifiableModifiedResidues = null;

        /**
         * Temporary save the amino acids for each call of identify().
         */
        private List<Group> residues;
	
	
	public ProteinModificationIdentifier(){
		
		bondLengthTolerance =  0.4;
		recordUnidentifiableModifiedCompounds = false;
		recordAdditionalAttachments = true;
		
		reset();
	}
	
	
	public void destroy(){
		if ( identifiedModifiedCompounds != null)
			identifiedModifiedCompounds.clear();
		if ( unidentifiableAtomLinkages != null)
			unidentifiableAtomLinkages.clear();
		if ( unidentifiableModifiedResidues != null)
			unidentifiableModifiedResidues.clear();
		
		unidentifiableAtomLinkages = null;
		unidentifiableAtomLinkages = null;
		unidentifiableModifiedResidues = null;
		
		
	}
	
	/**
	 * 
	 * @param bondLengthTolerance tolerance of error (in Angstroms) of the
	 *  covalent bond length, when calculating the atom distance threshold.
	 */
	public void setbondLengthTolerance(final double bondLengthTolerance) {
		if (bondLengthTolerance<0) {
			throw new IllegalArgumentException("bondLengthTolerance " +
					"must be positive.");
		}
		this.bondLengthTolerance = bondLengthTolerance;
	}
	
	/**
	 * 
	 * @param recordUnidentifiableAtomLinkages true if choosing to record unidentifiable
	 *  atoms; false, otherwise.
	 * @see #getRecordUnidentifiableCompounds
	 * @see #getUnidentifiableModifiedResidues
	 * @see #getUnidentifiableAtomLinkages
	 */
	public void setRecordUnidentifiableCompounds(boolean recordUnidentifiableModifiedCompounds) {
		this.recordUnidentifiableModifiedCompounds = recordUnidentifiableModifiedCompounds;
	}
	
	/**
	 * 
	 * @return true if choosing to record unidentifiable
	 *  atoms; false, otherwise.
	 * @see #setRecordUnidentifiableCompounds
	 * @see #getUnidentifiableModifiedResidues
	 * @see #getUnidentifiableAtomLinkages
	 */
	public boolean getRecordUnidentifiableCompounds() {
		return recordUnidentifiableModifiedCompounds;
	}
	
	/**
	 * 
	 * @param recordAdditionalAttachments true if choosing to record additional attachments
	 *  that are not directly attached to a modified residue.
	 * @see #getRecordAdditionalAttachments
	 */
	public void setRecordAdditionalAttachments(boolean recordAdditionalAttachments) {
		this.recordAdditionalAttachments = recordAdditionalAttachments;
	}
	
	/**
	 * 
	 * @return true if choosing to record additional attachments
	 *  that are not directly attached to a modified residue.
	 * @see #setRecordAdditionalAttachments
	 */
	public boolean getRecordAdditionalAttachments() {
		return recordAdditionalAttachments;
	}
	
	/**
	 * 
	 * @return a set of identified {@link ModifiedCompound}s from
	 *  the last parse result.
	 * @see ModifiedCompound
	 */
	public Set<ModifiedCompound> getIdentifiedModifiedCompound() {
		if (identifiedModifiedCompounds==null) {
			throw new IllegalStateException("No result available. Please call parse() first.");
		}
		
		return identifiedModifiedCompounds;
	}
	
	/**
	 * 
	 * @return a set of atom linkages, which represent the 
	 *  atom bonds that were not covered by the identified 
	 *  {@link ModifiedCompound}s from the last parse result.
	 *  Each element of the list is a array containing two atoms.
	 * @see StructureAtomLinkage
	 * @see #setRecordUnidentifiableCompounds
	 */
	public Set<StructureAtomLinkage> getUnidentifiableAtomLinkages() {
		if (!recordUnidentifiableModifiedCompounds) {
			throw new UnsupportedOperationException("Recording unidentified atom linkages" +
					"is not supported. Please setRecordUnidentifiableCompounds(true) first.");
		}
		
		if (identifiedModifiedCompounds==null) {
			throw new IllegalStateException("No result available. Please call parse() first.");
		}
		
		return unidentifiableAtomLinkages;
	}
	
	/**
	 * 
	 * @return a set of modified residues that were not covered by
	 *  the identified ModifiedCompounds from the last parse 
	 *  result.
	 *  @see StructureGroup
	 *  @see #setRecordUnidentifiableCompounds
	 *  @see #getIdentifiedModifiedCompound
	 */
	public Set<StructureGroup> getUnidentifiableModifiedResidues() {
		if (!recordUnidentifiableModifiedCompounds) {
			throw new UnsupportedOperationException("Recording unidentified atom linkages" +
					"is not supported. Please setRecordUnidentifiableCompounds(true) first.");
		}
		
		if (identifiedModifiedCompounds==null) {
			throw new IllegalStateException("No result available. Please call parse() first.");
		}
		
		return unidentifiableModifiedResidues;
	}
	
	/**
	 * Identify all registered modifications in a structure.
	 * @param structure
	 */
	public void identify(final Structure structure) {
		identify(structure, ProteinModificationRegistry.allModifications());
	}

	/**
	 * Identify a set of modifications in a structure.
	 * @param structure query {@link Structure}.
	 * @param potentialModifications query {@link ProteinModification}s.
	 */
	public void identify(final Structure structure,
			final Set<ProteinModification> potentialModifications) {
		if (structure==null) {
			throw new IllegalArgumentException("Null structure.");
		}
		
		identify(structure.getChains(), potentialModifications);
	}
	
	/**
	 * Identify all registered modifications in a chain. 
	 * @param chain query {@link Chain}.
	 */
	public void identify(final Chain chain) {
		identify(Collections.singletonList(chain));
	}
	
	/**
	 * Identify all registered modifications in chains. 
	 * @param chains query {@link Chain}s.
	 */
	public void identify(final List<Chain> chains) {
		identify(chains, ProteinModificationRegistry.allModifications());
	}
	
	/**
	 * Identify a set of modifications in a a chains.
	 * @param chain query {@link Chain}.
	 * @param potentialModifications query {@link ProteinModification}s.
	 */
	public void identify(final Chain chain,
			final Set<ProteinModification> potentialModifications)  {
		identify(Collections.singletonList(chain), potentialModifications);
	}
	
	/**
	 * Identify a set of modifications in a a list of chains.
	 * @param chains query {@link Chain}s.
	 * @param potentialModifications query {@link ProteinModification}s.
	 */
	public void identify(final List<Chain> chains,
			final Set<ProteinModification> potentialModifications) {
		
		if (chains==null) {
			throw new IllegalArgumentException("Null structure.");
		}
		
		if (potentialModifications==null) {
			throw new IllegalArgumentException("Null potentialModifications.");
		}
		
			
		reset();
		
		if (potentialModifications.isEmpty()) {
			return;
		}
		
		Map<String, Chain> mapChainIdChain = new HashMap<String, Chain>(chains.size());
		residues = new ArrayList<Group>();
		List<Group> ligands = new ArrayList<Group>();
		Map<Component, Set<Group>> mapCompGroups = new HashMap<Component, Set<Group>>();
		
		for (Chain chain : chains) {
			mapChainIdChain.put(chain.getChainID(), chain);
					
			List<Group> ress = StructureUtil.getAminoAcids(chain);
			List<Group> ligs = chain.getAtomLigands();
			residues.addAll(ress);
                        residues.removeAll(ligs);
			ligands.addAll(ligs);
			addModificationGroups(potentialModifications, ress, ligs, mapCompGroups);
		}
		
		if (residues.isEmpty()) {
			String pdbId = "?";
			if ( chains.size() > 0) {
				Structure struc = chains.get(0).getParent();
				if ( struc != null)
					pdbId = struc.getPDBCode(); 
			}
			System.err.println("WARNING: no amino acids found for "+ pdbId + ". Either you did not parse the PDB file with alignSEQRES records, or this record does not contain any amino acids.");
		}
		List<ModifiedCompound> modComps = new ArrayList<ModifiedCompound>();
		
		for (ProteinModification mod : potentialModifications) {
			ModificationCondition condition = mod.getCondition();
			List<Component> components = condition.getComponents();
			if (!mapCompGroups.keySet().containsAll(components)) {
				// not all components exist for this mod.
				continue;
			}
			
			int sizeComps = components.size();
			if (sizeComps==1) {
				
				processCrosslink1(mapCompGroups, modComps, mod, components);
			
			} else {
				
				processMultiCrosslink(mapCompGroups, modComps, mod, condition);
			}
		}

		if (recordAdditionalAttachments) {
			// identify additional groups that are not directly attached to amino acids.
			for (ModifiedCompound mc : modComps) {
				identifyAdditionalAttachments(mc, ligands, mapChainIdChain);
			}
		}
		
		mergeModComps(modComps);
		
		identifiedModifiedCompounds.addAll(modComps);
		
		
		// record unidentifiable linkage
		if (recordUnidentifiableModifiedCompounds) {
			recordUnidentifiableAtomLinkages(modComps, ligands);
			recordUnidentifiableModifiedResidues(modComps);
		}
	}

	private void reset() {
		identifiedModifiedCompounds = new LinkedHashSet<ModifiedCompound>();
		if (recordUnidentifiableModifiedCompounds) {
			unidentifiableAtomLinkages = new LinkedHashSet<StructureAtomLinkage>();
			unidentifiableModifiedResidues = new LinkedHashSet<StructureGroup>();
		}
		
	}

	private void processMultiCrosslink(
			Map<Component, Set<Group>> mapCompGroups,
			List<ModifiedCompound> modComps, ProteinModification mod,
			ModificationCondition condition) {
		// for multiple components
		
		// find linkages first
		List<List<Atom[]>> matchedAtomsOfLinkages =
				getMatchedAtomsOfLinkages(condition, mapCompGroups);
		
		if (matchedAtomsOfLinkages.size() != condition.getLinkages().size()) {
			return;
		} 
		
		assembleLinkages(matchedAtomsOfLinkages, mod, modComps);
		
	}

	private void processCrosslink1(Map<Component, Set<Group>> mapCompGroups,
			List<ModifiedCompound> modComps, ProteinModification mod,
			List<Component> components) {
		// modified residue
		// TODO: is this the correct logic for CROSS_LINK_1?
		Set<Group> modifiedResidues = mapCompGroups.get(components.get(0));
		if (modifiedResidues != null) {
			for (Group residue : modifiedResidues) {
				StructureGroup strucGroup = StructureUtil.getStructureGroup(residue, true);
				ModifiedCompound modRes = new ModifiedCompoundImpl(mod, strucGroup);
				modComps.add(modRes);
			}
		}
	}
	
	/**
	 * identify additional groups that are not directly attached to amino acids.
	 * @param mc {@link ModifiedCompound}.
	 * @param chain a {@link Chain}.
	 * @return a list of added groups.
	 */
	private void identifyAdditionalAttachments(ModifiedCompound mc, 
			List<Group> ligands, Map<String, Chain> mapChainIdChain) {
		if (ligands.isEmpty()) {
			return;
		}
		
		// TODO: should the additional groups only be allowed to the identified 
		// ligands or both amino acids and ligands? Currently only on ligands
		// ligands to amino acid bonds for same modification of unknown category
		// will be combined in mergeModComps()
		// TODO: how about chain-chain links?
		List<Group> identifiedGroups = new ArrayList<Group>();
		for (StructureGroup num : mc.getGroups(false)) {
			Group group;
			try {
				//String numIns = "" + num.getResidueNumber();
				//if (num.getInsCode() != null) {
				//	numIns += num.getInsCode();
				//}
				ResidueNumber resNum = new ResidueNumber();
				resNum.setChainId(num.getChainId());
				resNum.setSeqNum(num.getResidueNumber());
				resNum.setInsCode(num.getInsCode());
				//group = chain.getGroupByPDB(numIns);
				group = mapChainIdChain.get(num.getChainId()).getGroupByPDB(resNum);
			} catch (StructureException e) {
				e.printStackTrace();
				// should not happen
				continue;
			}
			identifiedGroups.add(group);
		}
		
		int start = 0;
		
		int n = identifiedGroups.size();
		while (n > start) {
			for (Group group1 : ligands) {
				for (int i=start; i<n; i++) {
					Group group2 = identifiedGroups.get(i);
					if (!identifiedGroups.contains(group1)) {
						List<Atom[]> linkedAtoms = StructureUtil.findAtomLinkages(
								group1, group2, false, bondLengthTolerance);
						if (!linkedAtoms.isEmpty()) {
							for (Atom[] atoms : linkedAtoms) {
								mc.addAtomLinkage(StructureUtil.getStructureAtomLinkage(atoms[0], 
										false, atoms[1], false));
							}
							identifiedGroups.add(group1);
							break;
						}
					}
				}
			}
			
			start = n;
			n = identifiedGroups.size();
		}
	}
	
	/**
	 * Merge identified modified compounds if linked.
	 */
	private void mergeModComps(List<ModifiedCompound> modComps) {
		TreeSet<Integer> remove = new TreeSet<Integer>();
		int n = modComps.size();
		for (int icurr=1; icurr<n; icurr++) {
			ModifiedCompound curr = modComps.get(icurr);
			
			String id = curr.getModification().getId();
			if (ProteinModificationRegistry.getById(id).getCategory()
					!=ModificationCategory.UNDEFINED)
				continue;
			
			// find linked compounds that before curr
			//List<Integer> merging = new ArrayList<Integer>();
			int ipre = 0;
			for (; ipre<icurr; ipre++) {
				if (remove.contains(ipre))	continue;
				ModifiedCompound pre = modComps.get(ipre);
				if (!Collections.disjoint(pre.getGroups(false),
						curr.getGroups(false))) {
					break;
				}
			}
			
			if (ipre<icurr) {				
				ModifiedCompound mcKeep = modComps.get(ipre);
				
				// merge modifications of the same type
				if (mcKeep.getModification().getId().equals(id)) {
					// merging the current one to the previous one
					mcKeep.addAtomLinkages(curr.getAtomLinkages());
					remove.add(icurr);
				}
			}
		}
		
		Iterator<Integer> it = remove.descendingIterator();
		while (it.hasNext()) {
			modComps.remove(it.next().intValue());
		}
	}
	
	/**
	 * Record unidentifiable atom linkages in a chain. Only linkages between two
	 * residues or one residue and one ligand will be recorded.
	 */
	private void recordUnidentifiableAtomLinkages(List<ModifiedCompound> modComps,
			List<Group> ligands) {
		
		// first put identified linkages in a map for fast query
		Set<StructureAtomLinkage> identifiedLinkages = new HashSet<StructureAtomLinkage>();
		for (ModifiedCompound mc : modComps) {
			identifiedLinkages.addAll(mc.getAtomLinkages());
		}
		
		// record
		// cross link
		int nRes = residues.size();
		for (int i=0; i<nRes-1; i++) {
			Group group1 = residues.get(i);
			for (int j=i+1; j<nRes; j++) {
				Group group2 = residues.get(j);
				List<Atom[]> linkages = StructureUtil.findAtomLinkages(
						group1, group2, true, bondLengthTolerance);
				for (Atom[] atoms : linkages) {
					StructureAtomLinkage link = StructureUtil.getStructureAtomLinkage(atoms[0], 
							true, atoms[1], true);
					unidentifiableAtomLinkages.add(link);
				}
			}
		}
		
		// attachment
		int nLig = ligands.size();
		for (int i=0; i<nRes; i++) {
			Group group1 = residues.get(i);
			for (int j=0; j<nLig; j++) {
				Group group2 = ligands.get(j);
				if (group1.equals(group2)) { // overlap between residues and ligands
					continue;
				}
				List<Atom[]> linkages = StructureUtil.findAtomLinkages(
						group1, group2, false, bondLengthTolerance);
				for (Atom[] atoms : linkages) {
					StructureAtomLinkage link = StructureUtil.getStructureAtomLinkage(atoms[0], 
							true, atoms[1], false);
					unidentifiableAtomLinkages.add(link);
				}
			}
		}
	}
	
	private void recordUnidentifiableModifiedResidues(List<ModifiedCompound> modComps) {
		Set<StructureGroup> identifiedComps = new HashSet<StructureGroup>();
		for (ModifiedCompound mc : modComps) {
			identifiedComps.addAll(mc.getGroups(true));
		}
		
		// TODO: use the ModifiedAminoAcid after Andreas add that.
		for (Group group : residues) {
			if (group.getType().equals(GroupType.HETATM)) {
				StructureGroup strucGroup = StructureUtil.getStructureGroup(
						group, true);
				if (!identifiedComps.contains(strucGroup)) {
					unidentifiableModifiedResidues.add(strucGroup);
				}
			}
		}
	}
	
	/**
	 * 
	 * @param modifications a set of {@link ProteinModification}s.
	 * @param residues
	 * @param ligands 
	 * @param saveTo save result to
	 * @return map from component to list of corresponding residues
	 *  in the chain.
	 */
	private void addModificationGroups(
			final Set<ProteinModification> modifications,
			final List<Group> residues,
			final List<Group> ligands,
			final Map<Component, Set<Group>> saveTo) {
		if (residues==null || ligands==null || modifications==null) {
			throw new IllegalArgumentException("Null argument(s).");
		}
		
		Map<Component,Set<Component>> mapSingleMultiComps = new HashMap<Component,Set<Component>>();
		for (ProteinModification mod : modifications) {
			ModificationCondition condition = mod.getCondition();
			for (Component comp : condition.getComponents()) {
				for (String pdbccId : comp.getPdbccIds()) {
					Component single = Component.of(Collections.singleton(pdbccId), 
							comp.isNTerminal(), comp.isCTerminal());
					Set<Component> mult = mapSingleMultiComps.get(single);
					if (mult == null) {
						mult = new HashSet<Component>();
						mapSingleMultiComps.put(single, mult);
					}
					mult.add(comp);
				}
			}
		}
		
		{
			// ligands
			Set<Component> ligandsWildCard = mapSingleMultiComps.get(
					Component.of("*"));
			for (Group group : ligands) {
				String pdbccId = group.getPDBName().trim();
				Set<Component> comps = mapSingleMultiComps.get(
						Component.of(pdbccId));
				
				for (Component comp : unionComponentSet(ligandsWildCard, comps)) {
					Set<Group> gs = saveTo.get(comp);
					if (gs==null) {
						gs = new LinkedHashSet<Group>();
						saveTo.put(comp, gs);
					}
					gs.add(group);
				}
			}
		}
		
		{
			// residues
			if (residues.isEmpty()) {
				return;
			}
			
			Set<Component> residuesWildCard = mapSingleMultiComps.get(
					Component.of("*"));
			
			// for all residues
			for (Group group : residues) {
				String pdbccId = group.getPDBName().trim();
				Set<Component> comps = mapSingleMultiComps.get(
						Component.of(pdbccId));
				
				for (Component comp : unionComponentSet(residuesWildCard, comps)) {
					Set<Group> gs = saveTo.get(comp);
					if (gs==null) {
						gs = new LinkedHashSet<Group>();
						saveTo.put(comp, gs);
					}
					gs.add(group);
				}
			}

			// for N-terminal
			int nRes = residues.size();
			int iRes = 0;
			Group res;
			do {
				// for all ligands on N terminal and the first residue
				res = residues.get(iRes++);

				Set<Component> nTermWildCard = mapSingleMultiComps.get(
						Component.of("*", true, false));

				Set<Component> comps = mapSingleMultiComps.get(
						Component.of(res.getPDBName(), true, false));
				
				for (Component comp : unionComponentSet(nTermWildCard, comps)) {
					Set<Group> gs = saveTo.get(comp);
					if (gs==null) {
						gs = new LinkedHashSet<Group>();
						saveTo.put(comp, gs);
					}
					gs.add(res);
				}
			} while (iRes<nRes && ligands.contains(res));
			
			// for C-terminal
			iRes = residues.size()-1;
			do {
				// for all ligands on C terminal and the last residue
				res = residues.get(iRes--);

				Set<Component> cTermWildCard = mapSingleMultiComps.get(
						Component.of("*", false, true));
				
				Set<Component> comps = mapSingleMultiComps.get(
						Component.of(res.getPDBName(), false, true));

				for (Component comp : unionComponentSet(cTermWildCard, comps)) {
					Set<Group> gs = saveTo.get(comp);
					if (gs==null) {
						gs = new LinkedHashSet<Group>();
						saveTo.put(comp, gs);
					}
					gs.add(res);
				}
			} while (iRes>=0 && ligands.contains(res));
		}
	}
	
	private Set<Component> unionComponentSet(Set<Component> set1, Set<Component> set2) {
		if (set1 == null && set2 == null)
			return Collections.emptySet();
		
		if (set1 == null)
			return set2;
		
		if (set2 == null)
			return set1;
		
		Set<Component> set = new HashSet<Component>(set1.size()+set2.size());
		set.addAll(set1);
		set.addAll(set2);
		
		return set;
	}
	
	/**
	 * Get matched atoms for all linkages.	
	 */
	private List<List<Atom[]>> getMatchedAtomsOfLinkages(
			ModificationCondition condition, Map<Component, Set<Group>> mapCompGroups) {
		List<ModificationLinkage> linkages = condition.getLinkages();
		int nLink = linkages.size();

		List<List<Atom[]>> matchedAtomsOfLinkages = 
				new ArrayList<List<Atom[]>>(nLink);
		
		for (int iLink=0; iLink<nLink; iLink++) {
			ModificationLinkage linkage = linkages.get(iLink);
			Component comp1 = linkage.getComponent1();
			Component comp2 = linkage.getComponent2();

//			boolean isAA1 = comp1.;
//			boolean isAA2 = comp2.getType()==true;
			
			Set<Group> groups1 = mapCompGroups.get(comp1);
			Set<Group> groups2 = mapCompGroups.get(comp2);
			
			List<Atom[]> list = new ArrayList<Atom[]>();

			List<String> potentialNamesOfAtomOnGroup1 = linkage.getPDBNameOfPotentialAtomsOnComponent1();
			for (String name : potentialNamesOfAtomOnGroup1) {
				if (name.equals("*")) {
					// wildcard
					potentialNamesOfAtomOnGroup1 = null; // search all atoms
					break;
				}
			}
			
			List<String> potentialNamesOfAtomOnGroup2 = linkage.getPDBNameOfPotentialAtomsOnComponent2();
			for (String name : potentialNamesOfAtomOnGroup2) {
				if (name.equals("*")) {
					// wildcard
					potentialNamesOfAtomOnGroup2 = null; // search all atoms
					break;
				}
			}

			for (Group g1 : groups1) {
				for (Group g2 : groups2) {
					if (g1.equals(g2)) {
						continue;
					}

                                        // only for wildcard match of two residues
                                        boolean ignoreNCLinkage = 
                                                potentialNamesOfAtomOnGroup1 == null &&
                                                potentialNamesOfAtomOnGroup2 == null &&
                                                residues.contains(g1) &&
                                                residues.contains(g2);
		
					Atom[] atoms = StructureUtil.findNearestAtomLinkage(
							g1, g2, 
							potentialNamesOfAtomOnGroup1,
							potentialNamesOfAtomOnGroup2,
                                                        ignoreNCLinkage,
							bondLengthTolerance);
					if (atoms!=null) {
						list.add(atoms);
					}
				}
			}
				
			if (list.isEmpty()) {
				// broken linkage
				break;
			}
	
			matchedAtomsOfLinkages.add(list);
		}
		
		return matchedAtomsOfLinkages;
	}
	
	/**
	 * Assembly the matched linkages.
	 * @param matchedAtomsOfLinkages
	 * @param mod
	 * @param condition
	 * @param ret ModifiedCompound will be stored here.
	 */
	private void assembleLinkages(List<List<Atom[]>> matchedAtomsOfLinkages,
			ProteinModification mod, List<ModifiedCompound> ret) {
		ModificationCondition condition = mod.getCondition();
		List<ModificationLinkage> modLinks = condition.getLinkages();
		
		int nLink = matchedAtomsOfLinkages.size();
		int[] indices = new int[nLink];
		Set<ModifiedCompound> identifiedCompounds = new HashSet<ModifiedCompound>();
		while (indices[0]<matchedAtomsOfLinkages.get(0).size()) {
			List<Atom[]> atomLinkages = new ArrayList<Atom[]>(nLink);
			for (int iLink=0; iLink<nLink; iLink++) {
				Atom[] atoms = matchedAtomsOfLinkages.get(iLink).get(indices[iLink]);
				atomLinkages.add(atoms);
			}
			if (matchLinkages(modLinks, atomLinkages)) {
				// matched
				
				int n = atomLinkages.size();
				List<StructureAtomLinkage> linkages = new ArrayList<StructureAtomLinkage>(n);
				for (int i=0; i<n; i++) {
					Atom[] linkage = atomLinkages.get(i);
					StructureAtomLinkage link = StructureUtil.getStructureAtomLinkage(
							linkage[0], residues.contains(linkage[0].getGroup()),
							linkage[1], residues.contains(linkage[1].getGroup()));
					linkages.add(link);
				}
				
				ModifiedCompound mc = new ModifiedCompoundImpl(mod, linkages);
				if (!identifiedCompounds.contains(mc)) {
					ret.add(mc);
					identifiedCompounds.add(mc);
				}
			}
			
			// indices++ (e.g. [0,0,1]=>[0,0,2]=>[1,2,0])
			int i = nLink-1;
			while (i>=0) {
				if (i==0 || indices[i]<matchedAtomsOfLinkages.get(i).size()-1) {
					indices[i]++;
					break;
				} else {
					indices[i] = 0;
					i--;
				}
			}
		}
	}
	
	/**
	 * 
	 * @param linkages
	 * @param atomLinkages
	 * @return true if atomLinkages satisfy the condition; false, otherwise.
	 */
	private boolean matchLinkages(List<ModificationLinkage> linkages, 
			List<Atom[]> atomLinkages) {
		int nLink = linkages.size();
		if (nLink != atomLinkages.size()) {
			return false;
		}
		for (int i=0; i<nLink-1; i++) {
			ModificationLinkage link1 = linkages.get(i);
			Atom[] atoms1 = atomLinkages.get(i);
			for (int j=i+1; j<nLink; j++) {
				ModificationLinkage link2 = linkages.get(j);
				Atom[] atoms2 = atomLinkages.get(j);
				
				// check components
				if (((link1.getIndexOfComponent1()==link2.getIndexOfComponent1())
							!= (atoms1[0].getGroup().equals(atoms2[0].getGroup())))
					|| ((link1.getIndexOfComponent1()==link2.getIndexOfComponent2())
							!= (atoms1[0].getGroup().equals(atoms2[1].getGroup())))
					|| ((link1.getIndexOfComponent2()==link2.getIndexOfComponent1())
							!= (atoms1[1].getGroup().equals(atoms2[0].getGroup())))
					|| ((link1.getIndexOfComponent2()==link2.getIndexOfComponent2())
							!= (atoms1[1].getGroup().equals(atoms2[1].getGroup())))) {
					return false;
				}
				
				// check atoms
				String label11 = link1.getLabelOfAtomOnComponent1();
				String label12 = link1.getLabelOfAtomOnComponent2();
				String label21 = link2.getLabelOfAtomOnComponent1();
				String label22 = link2.getLabelOfAtomOnComponent2();
				if ((label11!=null && label21!=null && label11.equals(label21))
							!= (atoms1[0].equals(atoms2[0]))
					 || (label11!=null && label22!=null && label11.equals(label22))
							!= (atoms1[0].equals(atoms2[1]))
					 || (label12!=null && label21!=null && label12.equals(label21))
							!= (atoms1[1].equals(atoms2[0]))
					 || (label12!=null && label22!=null && label12.equals(label22))
							!= (atoms1[1].equals(atoms2[1]))) {
					return false;
				}
			}
		}
		
		return true;
	}
}
