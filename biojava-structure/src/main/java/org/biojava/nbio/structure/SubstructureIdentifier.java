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
 * Created on December 19, 2013
 * Author: Douglas Myers-Turnbull
 */

package org.biojava.nbio.structure;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.Grid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This is the canonical way to identify a part of a structure.
 *
 * <p>The current syntax allows the specification of a set of residues from
 * the first model of a structure. Future versions may be extended to represent
 * additional properties.
 *
 * <p>Identifiers should adhere to the following specification, although some
 * additional forms may be tolerated where unambiguous for backwards compatibility.
 * <pre>
 * 		name          := pdbID
 * 		               | pdbID '.' chainID
 * 		               | pdbID '.' range
 * 		range         := range (',' range)?
 * 		               | chainID
 * 		               | chainID '_' resNum '-' resNum
 * 		pdbID         := [0-9][a-zA-Z0-9]{3}
 * 		chainID       := [a-zA-Z0-9]+
 * 		resNum        := [-+]?[0-9]+[A-Za-z]?
 * </pre>
 * For example:
 * <pre>
 * 		1TIM                            #whole structure
 * 		1tim                            #same as above
 * 		4HHB.C                          #single chain
 * 		3AA0.A,B                        #two chains
 * 		4GCR.A_1-40                     #substructure
 *      3iek.A_17-28,A_56-294,A_320-377 #substructure of 3 disjoint parts
 * </pre>
 * More options may be added to the specification at a future time.

 * @author dmyersturnbull
 * @author Spencer Bliven
 */
public class SubstructureIdentifier implements StructureIdentifier {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(SubstructureIdentifier.class);
	

	private final String pdbId;
	private final List<ResidueRange> ranges;

	/**
	 * Create a new identifier from a string.
	 * @param id
	 */
	public SubstructureIdentifier(String id) {
		String[] idRange = id.split("\\.");
		if(1 > idRange.length || idRange.length > 2 ) {
			throw new IllegalArgumentException(String.format("Malformed %s: %s",getClass().getSimpleName(),id));
		}
		if(idRange[0].length() != 4) {
			this.pdbId = idRange[0];
			// Changed from Exception to a warning to support files and stuff -sbliven 2015/01/22
			logger.warn(String.format("Unrecognized PDB code %s",this.pdbId));
		} else {
			this.pdbId = idRange[0].toUpperCase();
		}

		if( idRange.length == 2) {
			String rangeStr = idRange[1].trim();

			this.ranges = ResidueRange.parseMultiple(rangeStr);
		} else {
			this.ranges = new LinkedList<ResidueRange>();
		}
	}

	/**
	 * Create a new identifier based on a set of ranges.
	 *
	 * If ranges is empty, includes all residues.
	 * @param pdbId
	 * @param ranges
	 */
	public SubstructureIdentifier(String pdbId, List<ResidueRange> ranges) {
		if(ranges == null) {
			throw new NullPointerException("Null ranges list");
		}
		this.pdbId = pdbId;
		this.ranges = ranges;
	}

	@Override
	public String toString() {
		return getIdentifier();
	}

	/**
	 * Get the String form of this identifier.
	 *
	 * This provides the canonical form for a StructureIdentifier and has
	 * all the information needed to recreate a particular substructure.
	 *
	 * Example: 3iek.A_17-28,A_56-294
	 * @return The String form of this identifier
	 */
	@Override
	public String getIdentifier() {
		if (ranges.isEmpty()) return pdbId;
		return pdbId + "." + ResidueRange.toString(ranges);
	}

	public String getPdbId() {
		return pdbId;
	}

	public List<ResidueRange> getResidueRanges() {
		return ranges;
	}

	/**
	 * Return itself. SubstructureIdentifiers are canonical!
	 */
	@Override
	public SubstructureIdentifier toCanonical() {
		return this;
	}

	/**
	 * Takes a complete structure as input and reduces it to residues present in
	 * the specified ranges
	 *
	 * <p>The returned structure will be a shallow copy of the input, with shared
	 * Chains, Residues, etc.
	 * 
	 * <p>Ligands are handled in a special way. If a full chain is selected
	 * (e.g. '1ABC.A') then any waters and ligands with matching chain name are
	 * included. If a residue range is present ('1ABC.A:1-100') then any
	 * ligands (technically non-water non-polymer atoms) within
	 * {@link StructureTools#DEFAULT_LIGAND_PROXIMITY_CUTOFF} of the selected
	 * range are included, regardless of chain.
	 * @param input A full structure, e.g. as loaded from the PDB. The structure
	 * ID should match that returned by getPdbId().
	 * @return
	 * @throws StructureException
	 * @see StructureTools#getReducedStructure(Structure, String)
	 */
	@Override
	public Structure reduce(Structure s) throws StructureException {
		// Follows StructureImpl.clone()

		// Create new structure & copy basic properties
		Structure newS = new StructureImpl();

		newS.setPDBCode(s.getPDBCode());
		newS.setPDBHeader(s.getPDBHeader());
		newS.setName(this.toString());
		newS.setDBRefs(s.getDBRefs());
		newS.setBiologicalAssembly(s.isBiologicalAssembly());
		newS.getPDBHeader().setDescription(
				"sub-range " + ranges + " of "  + newS.getPDBCode() + " "
						+ s.getPDBHeader().getDescription());
		newS.setEntityInfos(new ArrayList<>());
		// TODO The following should be only copied for atoms which are present in the range.
		newS.setSSBonds(s.getSSBonds());
		newS.setSites(s.getSites());

		newS.setStructureIdentifier(this);

		for( int modelNr=0;modelNr<s.nrModels();modelNr++) {

			// Construct new model
			newS.addModel(new ArrayList<Chain>());

			if(getResidueRanges().isEmpty()) {
				// Include all residues
				newS.setEntityInfos(s.getEntityInfos());
				newS.setSSBonds(s.getSSBonds());
				newS.setSites(s.getSites());

				newS.setModel(modelNr, s.getModel(modelNr));
			} else {
				// Restrict residues
				for( ResidueRange range: getResidueRanges()) {

					String chainName = range.getChainName();
					ResidueNumber pdbresnum1 = range.getStart();
					ResidueNumber pdbresnum2 = range.getEnd();

//					StructureTools.addGroupsToStructure(newS, groups, modelNr, false);
					Chain polyChain; //polymer
					if(chainName.equals("_") ) {
						// Handle special case of "_" chain for single-chain proteins
						polyChain = s.getPolyChains(modelNr).get(0);
						chainName = polyChain.getName();
						if(pdbresnum1 != null)
							pdbresnum1.setChainName(chainName);
						if(pdbresnum2 != null)
							pdbresnum2.setChainName(chainName);

						if(s.getPolyChains().size() != 1) {
							// SCOP 1.71 uses this for some proteins with multiple chains
							// Print a warning in this ambiguous case
							logger.warn("Multiple possible chains match '_'. Using chain {}",chainName);
						}
					} else {
						// Explicit chain
						polyChain = s.getPolyChainByPDB(chainName,modelNr);
						if( polyChain == null ) {
							// Chain not found
							// Maybe it was a chain index, masquerading as a chainName?
							try {
								int chainNum = Integer.parseInt(chainName);
								polyChain = s.getChainByIndex(modelNr, chainNum);
								chainName = polyChain.getName();
								if(pdbresnum1 != null)
									pdbresnum1.setChainName(chainName);
								if(pdbresnum2 != null)
									pdbresnum2.setChainName(chainName);
								logger.warn("No chain found for {}. Interpretting it as an index, using chain {} instead",chainName,polyChain.getId());
							} catch(NumberFormatException e3) {
								// Not an index. Throw the original exception
								throw new StructureException(String.format("Unrecognized chain %s in %s",chainName,getIdentifier()));
							}
						}
					}

					if(pdbresnum1 == null && pdbresnum2 == null) {
						// Include all atoms with matching chainName
						StructureTools.addGroupsToStructure(newS, polyChain.getAtomGroups(), modelNr, false);
						for(Chain chain : s.getNonPolyChainsByPDB(chainName, modelNr) ) {
							StructureTools.addGroupsToStructure(newS, chain.getAtomGroups(), modelNr, false);
						}
						Chain waters = s.getWaterChainByPDB(chainName, modelNr);
						if( waters != null) {
							StructureTools.addGroupsToStructure(newS, waters.getAtomGroups(), modelNr, false);
						}
						
						// TODO do we need to prune SeqRes down to the atoms present? -SB 2016-10-7
					} else {
						// Include polymer range and any proximal ligands
						List<Group> polygroups = Arrays.asList(polyChain.getGroupsByPDB(pdbresnum1, pdbresnum2));
						StructureTools.addGroupsToStructure(newS, polygroups, modelNr, false);
						copyLigandsByProximity(s,newS, StructureTools.DEFAULT_LIGAND_PROXIMITY_CUTOFF, modelNr, modelNr);
					}
				} // end range
			}
		} // end modelNr

		return newS;
	}

	/**
	 * Loads the complete structure based on {@link #getPdbId()}.
	 *
	 * @param AtomCache A source of structures
	 * @return A Structure containing at least the atoms identified by this,
	 *  or null if no PDB ID is set
	 * @throws StructureException For errors loading and parsing the structure
	 * @throws IOException Errors reading the structure from disk
	 */
	@Override
	public Structure loadStructure(AtomCache cache) throws IOException, StructureException {
		String pdb = getPdbId();
		if(pdb == null)
			return null;
		return cache.getStructureForPdbId(pdb);
	}

	/**
	 * Supplements the reduced structure with ligands from the full structure based on
	 * a distance cutoff. Ligand groups are moved (destructively) from full to reduced
	 * if they fall within the cutoff of any atom in the reduced structure.
	 * The {@link StructureTools#DEFAULT_LIGAND_PROXIMITY_CUTOFF default cutoff}
	 * is used.
	 * @param full Structure containing all ligands
	 * @param reduced Structure with a subset of the polymer groups from full
	 * @see StructureTools#getLigandsByProximity(java.util.Collection, Atom[], double)
	 */
	protected static void copyLigandsByProximity(Structure full, Structure reduced) {
		// Normal case where all models should be copied from full to reduced
		assert full.nrModels() >= reduced.nrModels();
		for(int model = 0; model< reduced.nrModels(); model++) {
			copyLigandsByProximity(full, reduced, StructureTools.DEFAULT_LIGAND_PROXIMITY_CUTOFF, model, model);
		}
	}
	/**
	 * Supplements the reduced structure with ligands from the full structure based on
	 * a distance cutoff. Ligand groups are moved (destructively) from full to reduced
	 * if they fall within the cutoff of any atom in the reduced structure.
	 * @param full Structure containing all ligands
	 * @param reduced Structure with a subset of the polymer groups from full
	 * @param cutoff Distance cutoff (Å)
	 * @param fromModel source model in full
	 * @param toModel destination model in reduced
	 * @see StructureTools#getLigandsByProximity(java.util.Collection, Atom[], double)
	 */
	protected static void copyLigandsByProximity(Structure full, Structure reduced, double cutoff, int fromModel, int toModel) {
		// Geometric hashing of the reduced structure
		Grid grid = new Grid(cutoff);
		Atom[] nonwaters = StructureTools.getAllNonHAtomArray(reduced,true,toModel);
		if( nonwaters.length < 1 )
			return;
		grid.addAtoms(nonwaters);

		full.getNonPolyChains(fromModel).stream() //potential ligand chains
		.flatMap((chain) -> chain.getAtomGroups().stream() ) // potential ligand groups
		.filter( (g) -> !g.isWater() ) // ignore waters
		.filter( (g) -> !g.isPolymeric() ) // already shouldn't be polymeric, but filter anyways
		.filter( (g) -> grid.hasAnyContact(Calc.atomsToPoints(g.getAtoms())) ) // must contact reduced
		.sequential() // Keeps ligands from the same chain together if possible
		.reduce((Chain)null, // reduction updates the chain guess
				(guess, g ) -> {
					boolean wasAdded;
					try {
						// Check that it's not in reduced already
						wasAdded = reduced.findGroup(g.getChainId(),
								g.getResidueNumber().toString(), toModel) != null;
					} catch (StructureException e) {
						// not found
						wasAdded = false;
					}
					if( !wasAdded ) {
						// Add the ligand to reduced
						// note this is not idempotent, but it is synchronized on reduced
						logger.info("Adding ligand group {} {} by proximity",g.getPDBName(), g.getResidueNumber().toPDB());
						return StructureTools.addGroupToStructure(reduced, g, toModel, guess, false);
					}
					return guess;
				},
				// update to the new guess
				(oldGuess, newGuess) -> newGuess );
	}

}
