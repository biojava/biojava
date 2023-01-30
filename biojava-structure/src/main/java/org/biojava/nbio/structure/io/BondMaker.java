/*
 *                  BioJava development code
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
 * Created on Mar. 6, 2014
 *
 */
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.chem.ChemComp;
import org.biojava.nbio.structure.chem.ChemCompBond;
import org.biojava.nbio.structure.chem.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.util.PDBTemporaryStorageUtils.LinkRecord;
import org.rcsb.cif.model.ValueKind;
import org.rcsb.cif.schema.mm.StructConn;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Adds polymer bonds for peptides and nucleotides based on distance cutoffs and
 * intra-group (residue) bonds based on data from the Chemical Component Dictionary
 * to the Structure object.
 *
 * TODO the current implementation adds bonds to the first model only. This
 * should be sufficient for homogeneous models, but here are a few inhomogeneous models
 * in the PDB. A better handling of models should be considered in the future.
 *
 * @author Peter Rose
 * @author Ulysse Carion
 *
 */
public class BondMaker {
	private static final Logger logger = LoggerFactory.getLogger(BondMaker.class);

	/**
	 * The types of bonds that are read from struct_conn (type specified in field conn_type_id)
	 */
	public static final Set<String> BOND_TYPES_TO_PARSE;
	static {
		BOND_TYPES_TO_PARSE = new HashSet<>();
		BOND_TYPES_TO_PARSE.add("disulf");
		BOND_TYPES_TO_PARSE.add("covale");
		BOND_TYPES_TO_PARSE.add("covale_base");
		BOND_TYPES_TO_PARSE.add("covale_phosphate");
		BOND_TYPES_TO_PARSE.add("covale_sugar");
		BOND_TYPES_TO_PARSE.add("modres");
	}


	/**
	 * Maximum peptide (C - N) bond length considered for bond formation
	 */
	private static final double MAX_PEPTIDE_BOND_LENGTH = 1.8;
	/**
	 * Maximum nucleotide (P - O3') bond length considered for bond formation
	 */
	private static final double MAX_NUCLEOTIDE_BOND_LENGTH = 2.1;

	private final Structure structure;
	private final FileParsingParameters params;

	public BondMaker(Structure structure, FileParsingParameters params) {
		this.structure = structure;
		this.params = params;
	}

	/**
	 * Creates bond objects and corresponding references in Atom objects:
	 * <li>
	 * peptide bonds: inferred from sequence and distances
	 * </li>
	 * <li>
	 * nucleotide bonds: inferred from sequence and distances
	 * </li>
	 * <li>
	 * intra-group (residue) bonds: read from the chemical component dictionary, via {@link org.biojava.nbio.structure.chem.ChemCompProvider}
	 * </li>
	 */
	public void makeBonds() {
		logger.debug("Going to start making bonds");
		formPeptideBonds();
		formNucleotideBonds();
		formIntraResidueBonds();
		trimBondLists();
	}

	private void formPeptideBonds() {
		for (int modelInd=0; modelInd<structure.nrModels(); modelInd++){
			for (Chain chain : structure.getChains(modelInd)) {
				List<Group> groups = chain.getSeqResGroups();

				for (int i = 0; i < groups.size() - 1; i++) {
					if (!(groups.get(i) instanceof AminoAcidImpl)
							|| !(groups.get(i + 1) instanceof AminoAcidImpl))
						continue;

					AminoAcidImpl tail = (AminoAcidImpl) groups.get(i);
					AminoAcidImpl head = (AminoAcidImpl) groups.get(i + 1);

					// atoms with no residue number don't have atom information
					if (tail.getResidueNumber() == null
							|| head.getResidueNumber() == null) {
						continue;
					}

					formBondAltlocAware(tail, "C", head, "N", MAX_PEPTIDE_BOND_LENGTH, 1);
				}
			}
		}
	}

	private void formNucleotideBonds() {
		for (int modelInd=0; modelInd<structure.nrModels(); modelInd++){
			for (Chain chain : structure.getChains(modelInd)) {
				List<Group> groups = chain.getSeqResGroups();

				for (int i = 0; i < groups.size() - 1; i++) {
					if (!(groups.get(i) instanceof NucleotideImpl)
							|| !(groups.get(i + 1) instanceof NucleotideImpl))
						continue;

					NucleotideImpl tail = (NucleotideImpl) groups.get(i);
					NucleotideImpl head = (NucleotideImpl) groups.get(i + 1);

					// atoms with no residue number don't have atom information
					if (tail.getResidueNumber() == null
							|| head.getResidueNumber() == null) {
						continue;
					}

					formBondAltlocAware(head, "P", tail, "O3'", MAX_NUCLEOTIDE_BOND_LENGTH, 1);
				}
			}
		}
	}

	private void formIntraResidueBonds() {
		for (int modelInd=0; modelInd<structure.nrModels(); modelInd++){
			for (Chain chain : structure.getChains(modelInd)) {
				List<Group> groups = chain.getAtomGroups();
				for (Group mainGroup : groups) {
					// atoms with no residue number don't have atom information
					if (mainGroup.getResidueNumber() == null) {
						continue;
					}
					// Now add support for altLocGroup
					List<Group> totList = new ArrayList<>();
					totList.add(mainGroup);
					totList.addAll(mainGroup.getAltLocs());

					// Now iterate through this list
					for(Group group : totList){

						ChemComp aminoChemComp = ChemCompGroupFactory.getChemComp(group.getPDBName());
						logger.debug("chemcomp for residue {}-{} has {} atoms and {} bonds",
								group.getPDBName(), group.getResidueNumber(), aminoChemComp.getAtoms().size(), aminoChemComp.getBonds().size());

						for (ChemCompBond chemCompBond : aminoChemComp.getBonds()) {
							// note we don't check distance to make this call not too expensive
							formBondAltlocAware(group, chemCompBond.getAtomId1(),
									group, chemCompBond.getAtomId2(), -1, chemCompBond.getNumericalBondOrder());
						}
					}
				}
			}

		}
	}

	/**
	 * Form bond between atoms of the given names and groups, respecting alt loc rules to form bonds:
	 * no bonds between differently named alt locs (that are not the default alt loc '.')
	 * and multiple bonds for default alt loc to named alt loc.
	 * @param g1 first group
	 * @param name1 name of atom in first group
	 * @param g2 second group
	 * @param name2 name of atom in second group
	 * @param maxAllowedLength max length, if atoms distance above this length no bond will be added. If negative no check on distance is performed.
	 * @param bondOrder the bond order to be set in the created bond(s)
	 */
	private void formBondAltlocAware(Group g1, String name1, Group g2, String name2, double maxAllowedLength, int bondOrder) {
		List<Atom> a1s = getAtoms(g1, name1);
		List<Atom> a2s = getAtoms(g2, name2);

		if (a1s.isEmpty() || a2s.isEmpty()) {
			// some structures may be incomplete and not store info
			// about all of their atoms
			return;
		}

		for (Atom a1:a1s) {
			for (Atom a2:a2s) {
				if (a1.getAltLoc() != null && a2.getAltLoc()!=null &&
						a1.getAltLoc()!=' ' && a2.getAltLoc()!=' ' &&
						a1.getAltLoc() != a2.getAltLoc()) {
					logger.debug("Skipping bond between atoms with differently named alt locs {} (altLoc '{}') -- {} (altLoc '{}')",
							a1.toString(), a1.getAltLoc(), a2.toString(), a2.getAltLoc());
					continue;
				}
				if (maxAllowedLength<0) {
					// negative maxAllowedLength means we don't check distance and always add bond
					logger.debug("Forming bond between atoms {}-{} and {}-{} with bond order {}",
							a1.getPDBserial(), a1.getName(), a2.getPDBserial(), a2.getName(), bondOrder);
					new BondImpl(a1, a2, bondOrder);
				} else {
					if (Calc.getDistance(a1, a2) < maxAllowedLength) {
						logger.debug("Forming bond between atoms {}-{} and {}-{} with bond order {}. Distance is below {}",
								a1.getPDBserial(), a1.getName(), a2.getPDBserial(), a2.getName(), bondOrder, maxAllowedLength);
						new BondImpl(a1, a2, bondOrder);
					} else {
						logger.debug("Not forming bond between atoms {}-{} and {}-{} with bond order {}, because distance is above {}",
								a1.getPDBserial(), a1.getName(), a2.getPDBserial(), a2.getName(), bondOrder, maxAllowedLength);
					}
				}
			}
		}
	}

	/**
	 * Get all atoms (including possible alt locs) in given group that are name with the given atom name
	 * @param g the group
	 * @param name the atom name
	 * @return list of all atoms, or empty list if no atoms with the name
	 */
	private List<Atom> getAtoms(Group g, String name) {
		List<Atom> atoms = new ArrayList<>();
		List<Group> groupsWithAltLocs = new ArrayList<>();
		groupsWithAltLocs.add(g);
		groupsWithAltLocs.addAll(g.getAltLocs());
		for (Group group : groupsWithAltLocs) {
			Atom a = group.getAtom(name);
			// Check for deuteration
			if (a==null && name.startsWith("H")) {
				a = group.getAtom(name.replaceFirst("H", "D"));
				// Check it is actually deuterated
				if (a!=null && !a.getElement().equals(Element.D)){
					a=null;
				}
			}
			if (a!=null)
				atoms.add(a);
		}
		return atoms;
	}

	private void trimBondLists() {
		for (int modelInd=0; modelInd<structure.nrModels(); modelInd++){
			for (Chain chain : structure.getChains(modelInd)) {
				for (Group group : chain.getAtomGroups()) {
					for (Atom atom : group.getAtoms()) {
						if (atom.getBonds()!=null && atom.getBonds().size() > 0) {
							((ArrayList<Bond>) atom.getBonds()).trimToSize();
						}
					}
				}
			}
		}
	}

	/**
	 * Creates disulfide bond objects and references in the corresponding Atoms objects, given
	 * a list of {@link SSBondImpl}s parsed from a PDB file.
	 * @param disulfideBonds
	 */
	public void formDisulfideBonds(List<SSBondImpl> disulfideBonds) {
		for (SSBondImpl disulfideBond : disulfideBonds) {
			formDisulfideBond(disulfideBond);
		}
	}

	private void formDisulfideBond(SSBondImpl disulfideBond) {
		try {
			// The PDB format uses author chain ids to reference chains. But one author chain id corresponds to multiple asym ids,
			// thus we need to grab all the possible asym ids (poly and nonpoly) and then try to find the atoms
			// See issue https://github.com/biojava/biojava/issues/929
			Chain polyChain1 = structure.getPolyChainByPDB(disulfideBond.getChainID1());
			Chain polyChain2 = structure.getPolyChainByPDB(disulfideBond.getChainID2());
			List<Chain> nonpolyChains1 = structure.getNonPolyChainsByPDB(disulfideBond.getChainID1());
			List<Chain> nonpolyChains2 = structure.getNonPolyChainsByPDB(disulfideBond.getChainID2());

			List<String> allChainIds1 = new ArrayList<>();
			List<String> allChainIds2 = new ArrayList<>();
			if (polyChain1!=null) allChainIds1.add(polyChain1.getId());
			if (polyChain2!=null) allChainIds2.add(polyChain2.getId());
			if (nonpolyChains1!=null) nonpolyChains1.forEach(npc -> allChainIds1.add(npc.getId()));
			if (nonpolyChains2!=null) nonpolyChains2.forEach(npc -> allChainIds2.add(npc.getId()));

			Map<Integer, Atom> a = getAtomFromRecordTryMultipleChainIds("SG", "", disulfideBond.getResnum1(), disulfideBond.getInsCode1(), allChainIds1);

			Map<Integer, Atom> b = getAtomFromRecordTryMultipleChainIds("SG", "", disulfideBond.getResnum2(), disulfideBond.getInsCode2(), allChainIds2);

			for(int i=0; i<structure.nrModels(); i++){
				if(a.containsKey(i) && b.containsKey(i)){
					// TODO determine what the actual bond order of this bond is; for
					// now, we're assuming they're single bonds
					if(!a.get(i).equals(b.get(i))){
						Bond ssbond =  new BondImpl(a.get(i), b.get(i), 1);
						structure.addSSBond(ssbond);
					}
				}
			}


		} catch (StructureException e) {
			// Note, in Calpha only mode the CYS SG's are not present.
			if (! params.isParseCAOnly()) {
				logger.warn("Could not find atoms specified in SSBOND record: {}",disulfideBond.toString());
			} else {
				logger.debug("Could not find atoms specified in SSBOND record while parsing in parseCAonly mode.");
			}
		}
	}

	/**
	 * Creates bond objects from a LinkRecord as parsed from a PDB file
	 * @param linkRecord the PDB-format LINK record
	 */
	public void formLinkRecordBond(LinkRecord linkRecord) {
		// only work with atoms that aren't alternate locations
		if (linkRecord.getAltLoc1().equals(" ")
				|| linkRecord.getAltLoc2().equals(" "))
			return;

		try {
			// The PDB format uses author chain ids to reference chains. But one author chain id corresponds to multiple asym ids,
			// thus we need to grab all the possible asym ids (poly and nonpoly) and then try to find the atoms
			// See issue https://github.com/biojava/biojava/issues/943
			Chain polyChain1 = structure.getPolyChainByPDB(linkRecord.getChainID1());
			Chain polyChain2 = structure.getPolyChainByPDB(linkRecord.getChainID2());
			List<Chain> nonpolyChains1 = structure.getNonPolyChainsByPDB(linkRecord.getChainID1());
			List<Chain> nonpolyChains2 = structure.getNonPolyChainsByPDB(linkRecord.getChainID2());
			Chain waterChain1 = structure.getWaterChainByPDB(linkRecord.getChainID1());
			Chain waterChain2 = structure.getWaterChainByPDB(linkRecord.getChainID2());

			List<String> allChainIds1 = new ArrayList<>();
			List<String> allChainIds2 = new ArrayList<>();
			if (polyChain1!=null) allChainIds1.add(polyChain1.getId());
			if (polyChain2!=null) allChainIds2.add(polyChain2.getId());
			if (nonpolyChains1!=null) nonpolyChains1.forEach(npc -> allChainIds1.add(npc.getId()));
			if (nonpolyChains2!=null) nonpolyChains2.forEach(npc -> allChainIds2.add(npc.getId()));
			if (waterChain1!=null && linkRecord.getResName1().equals("HOH")) allChainIds1.add(waterChain1.getId());
			if (waterChain2!=null && linkRecord.getResName2().equals("HOH")) allChainIds2.add(waterChain2.getId());

			Map<Integer, Atom> a = getAtomFromRecordTryMultipleChainIds(linkRecord.getName1(), linkRecord.getAltLoc1(), linkRecord.getResSeq1(), linkRecord.getiCode1(), allChainIds1);

			Map<Integer, Atom> b = getAtomFromRecordTryMultipleChainIds(linkRecord.getName2(), linkRecord.getAltLoc2(), linkRecord.getResSeq2(), linkRecord.getiCode2(), allChainIds2);

			for(int i=0; i<structure.nrModels(); i++){
				if(a.containsKey(i) && b.containsKey(i)){
					// TODO determine what the actual bond order of this bond is; for
					// now, we're assuming they're single bonds
					if(!a.get(i).equals(b.get(i))){
						new BondImpl(a.get(i), b.get(i), 1);
					}
				}
			}
		} catch (StructureException e) {
			// Note, in Calpha only mode the link atoms may not be present.
			if (! params.isParseCAOnly()) {
				logger.warn("Could not find atoms specified in LINK record: {}",linkRecord.toString());
			} else {
				logger.debug("Could not find atoms specified in LINK record while parsing in parseCAonly mode.");
			}

		}
	}

	private Map<Integer, Atom> getAtomFromRecordTryMultipleChainIds(String name, String altLoc, String resSeq, String iCode, List<String> chainIds) throws StructureException {
		Map<Integer, Atom> a = null;
		for (String chainId : chainIds) {
			try {
				a = getAtomFromRecord(name, altLoc, chainId, resSeq, iCode);
				// first instance that doesn't give an exception will be considered the right one. Not much more we can do here
				break;
			} catch (StructureException e) {
				logger.debug("Tried to get atom {} {} {} (alt loc {}) from chain id {}, but did not find it", name, resSeq, iCode, altLoc, chainId);
			}
		}
		if (a == null) {
			throw new StructureException("Could not find atom "+name+" "+resSeq+" "+iCode+" (alt loc "+altLoc+")");
		}
		return a;
	}


	public void formBondsFromStructConn(StructConn conn) {
		final String symop = "1_555"; // For now - accept bonds within origin asymmetric unit.
		List<Bond> ssbonds = new ArrayList<>();

		for (int i = 0; i < conn.getRowCount(); i++) {
			if (!BOND_TYPES_TO_PARSE.contains(conn.getConnTypeId().get(i))) continue;
			String chainId1;
			String chainId2;

			chainId1 = conn.getPtnr1LabelAsymId().get(i);
			chainId2 = conn.getPtnr2LabelAsymId().get(i);

			String insCode1 = "";
			if (conn.getPdbxPtnr1PDBInsCode().getValueKind(i) == ValueKind.PRESENT) {
				insCode1 = conn.getPdbxPtnr1PDBInsCode().get(i);
			}
			String insCode2 = "";
			if (conn.getPdbxPtnr2PDBInsCode().getValueKind(i) == ValueKind.PRESENT) {
				insCode2 = conn.getPdbxPtnr2PDBInsCode().get(i);
			}

			String seqId1 = conn.getPtnr1AuthSeqId().getStringData(i);
			String seqId2 = conn.getPtnr2AuthSeqId().getStringData(i);
			String resName1 = conn.getPtnr1LabelCompId().get(i);
			String resName2 = conn.getPtnr2LabelCompId().get(i);
			String atomName1 = conn.getPtnr1LabelAtomId().get(i);
			String atomName2 = conn.getPtnr2LabelAtomId().get(i);
			String altLoc1 = "";
			if (conn.getPdbxPtnr1LabelAltId().getValueKind(i) == ValueKind.PRESENT) {
				altLoc1 = conn.getPdbxPtnr1LabelAltId().get(i);
			}
			String altLoc2 = "";
			if (conn.getPdbxPtnr2LabelAltId().getValueKind(i) == ValueKind.PRESENT) {
				altLoc2 = conn.getPdbxPtnr2LabelAltId().get(i);
			}

			// TODO: when issue 220 is implemented, add robust symmetry handling to allow bonds between symmetry-related molecules.
			if (!conn.getPtnr1Symmetry().get(i).equals(symop) || !conn.getPtnr2Symmetry().get(i).equals(symop) ) {
				logger.info("Skipping bond between atoms {}(residue {}{}) and {}(residue {}{}) belonging to different symmetry partners, because it is not supported yet",
						atomName1, seqId1, insCode1, atomName2, seqId2, insCode2);
				continue;
			}


			String altLocStr1 = altLoc1.isEmpty()? "" : "(alt loc "+altLoc1+")";
			String altLocStr2 = altLoc2.isEmpty()? "" : "(alt loc "+altLoc2+")";

			Map<Integer,Atom> a1 = null;
			Map<Integer,Atom> a2 = null;

			try {
				a1 = getAtomFromRecord(atomName1, altLoc1, chainId1, seqId1, insCode1);

			} catch (StructureException e) {

				logger.warn("Could not find atom specified in struct_conn record: {}{}({}) in chain {}, atom {} {}", seqId1, insCode1, resName1, chainId1, atomName1, altLocStr1);
				continue;
			}
			try {
				a2 = getAtomFromRecord(atomName2, altLoc2, chainId2, seqId2, insCode2);
			} catch (StructureException e) {

				logger.warn("Could not find atom specified in struct_conn record: {}{}({}) in chain {}, atom {} {}", seqId2, insCode2, resName2, chainId2, atomName2, altLocStr2);
				continue;
			}

			if (a1==null) {
				// we couldn't find the atom, something must be wrong with the file
				logger.warn("Could not find atom {} {} from residue {}{}({}) in chain {} to create bond specified in struct_conn", atomName1, altLocStr1, seqId1, insCode1, resName1, chainId1);
				continue;
			}
			if (a2==null) {
				// we couldn't find the atom, something must be wrong with the file
				logger.warn("Could not find atom {} {} from residue {}{}({}) in chain {} to create bond specified in struct_conn", atomName2, altLocStr2, seqId2, insCode2, resName2, chainId2);
				continue;
			}

			// assuming order 1 for all bonds, no information is provided by struct_conn
			for (int j = 0; j < structure.nrModels(); j++) {
				Bond bond = null;
				if (a1.containsKey(j) && a2.containsKey(j)) {
					if (!a1.get(j).equals(a2.get(j))) {
						bond = new BondImpl(a1.get(j), a2.get(j), 1);
					}
				}
				if (bond != null) {
					if (conn.getConnTypeId().get(i).equals("disulf")) {
						ssbonds.add(bond);
					}
				}
			}
		}

		// only for ss bonds we add a specific map in structure, all the rests are linked only from Atom.getBonds
		structure.setSSBonds(ssbonds);
	}

	private Map<Integer,Atom> getAtomFromRecord(String name, String altLoc, String chainID, String resSeq, String iCode)
			throws StructureException {

		if (iCode==null || iCode.isEmpty()) {
			iCode = " "; // an insertion code of ' ' is ignored
		}
		Map<Integer, Atom> outMap = new HashMap<>();
		ResidueNumber resNum = new ResidueNumber(chainID, Integer.parseInt(resSeq), iCode.charAt(0));

		for (int i=0; i<structure.nrModels(); i++){
			Chain chain = structure.getChain(chainID,i);
			Group group = chain.getGroupByPDB(resNum);

			Group g = group;
			// there is an alternate location
			if (!altLoc.isEmpty()) {
				g = group.getAltLocGroup(altLoc.charAt(0));
				if (g==null)
					throw new StructureException("Could not find altLoc code "+altLoc+" in group "+resSeq+iCode+" of chain "+ chainID);
			}
			Atom a = g.getAtom(name);
			if (a!=null){
				outMap.put(i,a);
			}
		}
		return outMap;
	}
}
