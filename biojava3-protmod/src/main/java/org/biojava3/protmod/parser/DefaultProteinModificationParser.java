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

package org.biojava3.protmod.parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ComponentType;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModificationLinkage;
import org.biojava3.protmod.ModifiedCompound;
import org.biojava3.protmod.ModifiedCompoundImpl;
import org.biojava3.protmod.ProteinModification;

/**
 * Identify attachment modification in a 3-D structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class DefaultProteinModificationParser
implements ProteinModificationParser {
	
	private double bondLengthTolerance = 0.4;
	private boolean recordUnidentifiableAtomLinkages = false;
	
	private List<ModifiedCompound> identifiedModifiedCompounds = null;
	private List<Atom[]> unidentifiableAtomLinkages = null;
	
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
	 */
	public void setRecordUnidentifiableAtomLinkages(boolean recordUnidentifiableAtomLinkages) {
		this.recordUnidentifiableAtomLinkages = recordUnidentifiableAtomLinkages;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public List<ModifiedCompound> getIdentifiedModifiedCompound() {
		if (identifiedModifiedCompounds==null) {
			throw new IllegalStateException("No result available. Please call parse() first.");
		}
		
		return identifiedModifiedCompounds;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public List<Atom[]> getUnidentifiableAtomLinkages() {
		if (!recordUnidentifiableAtomLinkages) {
			throw new UnsupportedOperationException("Recording unidentified atom linkages" +
					"is not supported. Please setRecordUnidentifiableAtoms(true) first.");
		}
		
		if (identifiedModifiedCompounds==null) {
			throw new IllegalStateException("No result available. Please call parse() first.");
		}
		
		return unidentifiableAtomLinkages;
	}
	
	/**
	 * {@inheritDoc}
	 * @throws IllegalArgumentException if null structure or
	 *  potentialModification, or , or if the nodelnr is less then 0 or
	 *  larger than or equal to the number of models in the structure.
	 */
	@Override
	public void parse(final Structure structure, 
			final Set<ProteinModification> potentialModifications,
			final int modelnr) {
		identifiedModifiedCompounds = new ArrayList<ModifiedCompound>();
		if (recordUnidentifiableAtomLinkages) {
			unidentifiableAtomLinkages = new ArrayList<Atom[]>();
		}
		
		if (structure==null) {
			throw new IllegalArgumentException("Null structure.");
		}
		
		if (potentialModifications==null) {
			throw new IllegalArgumentException("Null potentialModifications.");
		}
		
		if (modelnr >= structure.nrModels() || modelnr < 0) {
			throw new IllegalArgumentException("modelnr should be between 0 to "
					+ structure.nrModels() + " for the structure "
					+ structure.getName());
		}
		
		if (potentialModifications.isEmpty()) {
			return;
		}
		
		List<Chain> chains = structure.getChains(modelnr);
		int nChains = chains.size();
		
		Map<Chain, List<Group>> mapChainResidues = new LinkedHashMap<Chain,List<Group>>(nChains);
		Map<Chain, List<Group>> mapChainLigands = new LinkedHashMap<Chain, List<Group>>(nChains);
		
		for (int i=0; i<nChains; i++) {
			Chain chain = chains.get(i);
			List<Group> residues = chain.getSeqResGroups();
			mapChainResidues.put(chain, residues);
			List<Group> ligands = new ArrayList<Group>(chain.getAtomLigands());
			if (ligands.removeAll(residues)) {
				System.err.println(structure.getPDBCode()+"\t"+modelnr+"\t"
						+chain.getName()+": Overlapping between ligands and residues");
			}
			mapChainLigands.put(chain, ligands);
		}
		
		for (Chain chain : chains) {
			Map<Component, List<Group>> mapCompGroups = 
				getModificationGroups(potentialModifications, 
						mapChainResidues.get(chain), mapChainLigands.get(chain));
			
			for (ProteinModification mod : potentialModifications) {
				ModificationCondition condition = mod.getCondition();
				List<Component> components = condition.getComponents();
				if (!mapCompGroups.keySet().containsAll(components)) {
					// not all components exist for this mod.
					continue;
				}
				
				int sizeComps = components.size();
				if (sizeComps==1) {
					// modified residue
					// TODO: is this the correct logic for CROSS_LINK_1?
					List<Group> residues = mapCompGroups.get(components.get(0));
					if (residues != null) {
						for (Group residue : residues) {
							ModifiedCompound modRes = new ModifiedCompoundImpl(mod, residue);
							identifiedModifiedCompounds.add(modRes);
						}
					}
				} else {
					// for multiple components
					
					// find linkages first
					List<List<Atom[]>> matchedAtomsOfLinkages =
							getMatchedAtomsOfLinkages(condition, mapCompGroups);
					
					if (matchedAtomsOfLinkages.size() != condition.getLinkages().size()) {
						continue;
					}
					
					assembleLinkages(matchedAtomsOfLinkages, mod, identifiedModifiedCompounds);
				}
			}
		}

		// identify additional groups that are not directly attached to amino acids. 
		for (ModifiedCompound mc : identifiedModifiedCompounds) {
			Chain chain = mc.getGroups().iterator().next().getParent();
			Set<Group> ligands = new LinkedHashSet<Group>(mapChainLigands.get(chain));
			identifyAdditionalAttachments(mc, ligands);
		}
		
		// record unidentifiable linkage
		if (recordUnidentifiableAtomLinkages) {
			for (Chain chain : chains) {
				List<Group> groups = new ArrayList<Group>();
				groups.addAll(mapChainResidues.get(chain));
				groups.addAll(mapChainLigands.get(chain));
				recordUnidentifiableAtomLinkages(chain, groups);
			}
		}
	}
	
	/**
	 * identify additional groups that are not directly attached to amino acids.
	 * @param mc {@link ModifiedCompound}.
	 * @return a list of added groups.
	 */
	private void identifyAdditionalAttachments(ModifiedCompound mc, 
			Set<Group> ligands) {
		if (ligands.isEmpty()) {
			return;
		}
		
		// TODO: should the additional groups only be allowed to the identified 
		// heta groups or both amino acids and heta groups?
		List<Group> identifiedLigands = new ArrayList<Group>();
		for (Group group : mc.getGroups()) {
			if (ligands.contains(group)) {
				identifiedLigands.add(group);
			}
		}
		
		int start = 0;
		
		int n = identifiedLigands.size();
		while (n > start) {
			for (Group group1 : ligands) {
				for (int i=start; i<n; i++) {
					Group group2 = identifiedLigands.get(i);
					if (!identifiedLigands.contains(group1)) {
						List<Atom[]> linkages = ProteinModificationParserUtil.
								findNonNCAtomLinkages(group2, group1, bondLengthTolerance);
						if (!linkages.isEmpty()) {
							for (Atom[] linkage : linkages) {
								mc.addAtomLinkage(linkage[0], linkage[1]);
							}
							identifiedLigands.add(group1);
							break;
						}
					}
				}
			}
			
			start = n;
			n = identifiedLigands.size();
		}
	}
	
	/**
	 * record unidentifiable atom linkages in a chain.
	 * @param chain a {@link Chain}.
	 */
	private void recordUnidentifiableAtomLinkages(Chain chain, List<Group> groups) {
		// first put identified linkages in a map for fast query
		Map<Atom, Set<Atom>> identifiedLinkages = new HashMap<Atom, Set<Atom>>();
		for (ModifiedCompound mc : identifiedModifiedCompounds) {
			List<Atom[]> linkages = mc.getAtomLinkages();
			for (Atom[] linkage : linkages) {
				Set<Atom> set = identifiedLinkages.get(linkage[0]);
				if (set == null) {
					set = new HashSet<Atom>();
					identifiedLinkages.put(linkage[0], set);
				}
				set.add(linkage[1]);
				
				set = identifiedLinkages.get(linkage[1]);
				if (set == null) {
					set = new HashSet<Atom>();
					identifiedLinkages.put(linkage[1], set);
				}
				set.add(linkage[0]);
			}
		}
		
		// record
		int n = groups.size();
		for (int i=0; i<n-1; i++) {
			Group group1 = groups.get(i);
			for (int j=i+1; j<n; j++) {
				Group group2 = groups.get(j);
				List<Atom[]> linkages = ProteinModificationParserUtil.
						findNonNCAtomLinkages(group1, group2, bondLengthTolerance);
				for (Atom[] linkage : linkages) {
					Set<Atom> set = identifiedLinkages.get(linkage[0]);
					if (set == null || !set.contains(linkage[1])) {
						unidentifiableAtomLinkages.add(linkage);
					}
				}
			}
		}
	}
	
	/**
	 * 
	 * @param modifications a set of {@link ProteinModification}s.
	 * @param residues
	 * @param ligands 
	 * @return map from component to list of corresponding residues
	 *  in the chain.
	 */
	private Map<Component, List<Group>> getModificationGroups(
			final Set<ProteinModification> modifications,
			final List<Group> residues,
			final List<Group> ligands) {
		if (residues==null || ligands==null || modifications==null) {
			throw new IllegalArgumentException("Null argument(s).");
		}
		
		Set<Component> comps = new HashSet<Component>();
		for (ProteinModification mod : modifications) {
			ModificationCondition condition = mod.getCondition();
			for (Component comp : condition.getComponents()) {
				comps.add(comp);
			}
		}
		
		Map<Component, List<Group>> mapCompRes = 
			new HashMap<Component, List<Group>>();
		
		{
			// ligands
			for (Group group : ligands) {
				String pdbccId = group.getPDBName().trim();
				Component comp = Component.of(pdbccId, ComponentType.LIGAND);
				if (!comps.contains(comp)) {
					continue;
				}
				List<Group> gs = mapCompRes.get(comp);
				if (gs==null) {
					gs = new ArrayList<Group>();
					mapCompRes.put(comp, gs);
				}
				gs.add(group);
			}
		}
		
		{
			// residues
			if (residues.isEmpty()) {
				return mapCompRes;
			}
			
			// for all residues
			for (Group group : residues) {
				String pdbccId = group.getPDBName().trim();
				Component comp = Component.of(pdbccId, ComponentType.AMINOACID);
				if (!comps.contains(comp)) {
					continue;
				}
				List<Group> gs = mapCompRes.get(comp);
				if (gs==null) {
					gs = new ArrayList<Group>();
					mapCompRes.put(comp, gs);
				}
				gs.add(group);
			}

			// for N-terminal
			Group res = residues.get(0);
			Component comp = Component.of(res.getPDBName(), ComponentType.AMINOACID, true, false);
			if (comps.contains(comp)) {
				List<Group> gs = Collections.singletonList(res);
				mapCompRes.put(comp, gs);
			}
			
			// for C-terminal
			res = residues.get(residues.size()-1);
			comp = Component.of(res.getPDBName(), ComponentType.AMINOACID, false, true);
			if (comps.contains(comp)) {
				List<Group> gs = Collections.singletonList(res);
				mapCompRes.put(comp, gs);
			}
		}

		return mapCompRes;
	}
	
	/**
	 * Get matched atoms for all linkages.	
	 */
	private List<List<Atom[]>> getMatchedAtomsOfLinkages(
			ModificationCondition condition, Map<Component, List<Group>> mapCompGroups) {
		List<ModificationLinkage> linkages = condition.getLinkages();
		int nLink = linkages.size();

		List<List<Atom[]>> matchedAtomsOfLinkages = 
				new ArrayList<List<Atom[]>>(nLink);
		
		for (int iLink=0; iLink<nLink; iLink++) {
			ModificationLinkage linkage = linkages.get(iLink);
			Component comp1 = linkage.getComponent1();
			Component comp2 = linkage.getComponent2();
			List<Group> groups1 = mapCompGroups.get(comp1);
			List<Group> groups2 = mapCompGroups.get(comp2);						
			
			List<Atom[]> list = new ArrayList<Atom[]>();

			List<String> potentialNamesOfAtomOnGroup1 = linkage.getPDBNameOfPotentialAtomsOnComponent1();
			List<String> potentialNamesOfAtomOnGroup2 = linkage.getPDBNameOfPotentialAtomsOnComponent2();

			for (Group g1 : groups1) {
				for (Group g2 : groups2) {
					if (g1 == g2) {
						continue;
					}
		
					Atom[] atoms = ProteinModificationParserUtil.findNearestNonNCAtomLinkage(g1, g2, 
							potentialNamesOfAtomOnGroup1, potentialNamesOfAtomOnGroup2, bondLengthTolerance);
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
		
		int nLink = matchedAtomsOfLinkages.size();
		int[] indices = new int[nLink];
		Set<ModifiedCompound> identifiedCompounds = new HashSet<ModifiedCompound>();
		while (indices[0]<matchedAtomsOfLinkages.get(0).size()) {
			List<Atom[]> atomLinkages = new ArrayList<Atom[]>(nLink);
			Set<Group> groups = new LinkedHashSet<Group>();
			for (int iLink=0; iLink<nLink; iLink++) {
				Atom[] atoms = matchedAtomsOfLinkages.get(iLink).get(indices[iLink]);
				atomLinkages.add(atoms);
				groups.add(atoms[0].getParent());
				groups.add(atoms[1].getParent());
			}
			if (matchLinkages(condition.getLinkages(), atomLinkages)) {
				// matched
				ModifiedCompound mc = new ModifiedCompoundImpl(mod, groups, atomLinkages);
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
	 * @param condition
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
							!= (atoms1[0].getParent()==atoms2[0].getParent()))
					|| ((link1.getIndexOfComponent1()==link2.getIndexOfComponent2())
							!= (atoms1[0].getParent()==atoms2[1].getParent()))
					|| ((link1.getIndexOfComponent2()==link2.getIndexOfComponent1())
							!= (atoms1[1].getParent()==atoms2[0].getParent()))
					|| ((link1.getIndexOfComponent2()==link2.getIndexOfComponent2())
							!= (atoms1[1].getParent()==atoms2[1].getParent()))) {
					return false;
				}
				
				// check atoms
				String label11 = link1.getLabelOfAtomOnComponent1();
				String label12 = link1.getLabelOfAtomOnComponent2();
				String label21 = link2.getLabelOfAtomOnComponent1();
				String label22 = link2.getLabelOfAtomOnComponent2();
				if ((label11!=null && label21!=null && label11.equals(label21))
							!= (atoms1[0]==atoms2[0])
					 || (label11!=null && label22!=null && label11.equals(label22))
							!= (atoms1[0]==atoms2[1])
					 || (label12!=null && label21!=null && label12.equals(label21))
							!= (atoms1[1]==atoms2[0])
					 || (label12!=null && label22!=null && label12.equals(label22))
							!= (atoms1[1]==atoms2[1])) {
					return false;
				}
			}
		}
		
		return true;
	}
}
