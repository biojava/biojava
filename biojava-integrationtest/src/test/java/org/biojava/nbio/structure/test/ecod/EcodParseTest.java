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
package org.biojava.nbio.structure.test.ecod;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomPositionMap;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.ResidueRangeAndLength;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.ecod.EcodDatabase;
import org.biojava.nbio.structure.ecod.EcodDomain;
import org.biojava.nbio.structure.ecod.EcodFactory;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This is not a unit test.
 *
 * It is a long-running parsing test, which sequentially parses all ECOD domains.
 *
 * The most common warning is caused by residue ranges with missing terminal CA atoms,
 * which cause a warning to print.
 *
 * develop134, 127, 84 and earlier versions also had a number of invalid ranges, which
 * cause error messages to print.
 *
 * Filtering log4j messages to the 'error' level will filter all but the most
 * grievous errors.
 *
 * Faster unit tests go in {@link EcodInstallationTest}.
 *
 * @author blivens
 *
 */
public class EcodParseTest {
	private static final Logger logger = LoggerFactory.getLogger(EcodParseTest.class);

	public static void main(String[] args) throws IOException {
//		String ecodVersion = "develop124";
		String ecodVersion = "latest";

		int errors = testVersion(ecodVersion);
		logger.info("Done. {} errors.",errors);

	}

	private static int testVersion(String ecodVersion) throws IOException {
		EcodDatabase ecod = EcodFactory.getEcodDatabase(ecodVersion);
		AtomCache cache = new AtomCache();
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		List<EcodDomain> domains = ecod.getAllDomains();
//		domains = Arrays.asList(ecod.getDomainsById("e1yfbB2"));
//		domains = Arrays.asList(ecod.getDomainsById("e1w50A2"));
//		domains = Arrays.asList(ecod.getDomainsById("e2ftlE1"));
		int errors = 0;
		for(EcodDomain d : domains) {
			Atom[] ca1;
			Structure struct;
			try {
				struct = cache.getStructure(d.getPdbId());
				ca1 = StructureTools.getRepresentativeAtomArray(struct);
			} catch (IOException e) {
				logger.error("Error getting structure for "+d.getDomainId(),e);
				errors++;
				continue;
			} catch (StructureException e) {
				logger.error("Error getting structure for "+d.getDomainId(),e);
				errors++;
				continue;
			}

			// Test that the ranges can be parsed
			String rangeStr = d.getRange();
			AtomPositionMap map = new AtomPositionMap(ca1);
			List<? extends ResidueRange> ranges;
			try {
				// Parses range given in domain
				ranges = ResidueRange.parseMultiple(rangeStr);
			} catch(Exception e) {
				logger.error("Error parsing "+d.getPdbId()+"_"+d.getRange(),e);
				errors++;
				continue;
			}
			boolean clean = true;
			for(ResidueRange r : ranges) {
				if( r == null ) {
					clean = false;
				}
			}
			if( ! clean ) {
				logger.error("Empty range for {}_{}",d.getPdbId(),d.getRange());
				errors++;
				continue;
			}


			// Check that the ranges are valid (or at least that they have a group)
			for(ResidueRange range : ranges) {
				try {
					Integer start = map.getPosition(range.getStart());
					if(start == null) {
						Group g = struct.getPolyChainByPDB(range.getStart().getChainName()).getGroupByPDB(range.getStart());
						if(g!=null) {
							logger.warn("No CA atom for starting residue "+d.getDomainId()+"_"+range);
							clean = false;
						} else {
							logger.error("Start doesn't exist for "+d.getDomainId()+"_"+range.toString());
							clean = false;
						}
					}
				} catch(Exception e) {
					logger.error("Start doesn't exist for "+d.getDomainId()+"_"+range.toString(),e);
					clean = false;
				}
				try {
					Integer end = map.getPosition(range.getEnd());
					if(end == null) {
						Group g = null;
						try {
							g = struct.getPolyChainByPDB(range.getEnd().getChainName()).getGroupByPDB(range.getEnd());
						} catch(StructureException e ) {}
						if(g!=null) {
							logger.warn("No CA atom for ending residue "+d.getDomainId()+"_"+range);
							clean = false;
						} else {
							logger.error("End doesn't exist for "+d.getDomainId()+"_"+range.toString());
							clean = false;
						}
					}
				} catch(Exception e) {
					logger.error("End doesn't exist for "+d.getDomainId()+"_"+range.toString(),e);
					clean = false;
				}
			}

			// Try to recover from missing residues by trimming them to the residue range
			try {
				// Parses more flexibly, giving only a warning for missing residues
				ranges = ResidueRangeAndLength.parseMultiple(rangeStr,map);
			} catch(Exception e) {
				logger.error("Error parsing "+d.getPdbId()+"_"+d.getRange(),e);
				errors++;
				continue;
			}
			clean = true;
			for(ResidueRange r : ranges) {
				if( r == null ) {
					clean = false;
				}
			}
			if( ! clean ) {
				logger.error("Empty range for {}_{}",d.getPdbId(),d.getRange());
				errors++;
				continue;
			}

			// Test whether we can use it to get a structure
			String pdbRangeStr = String.format("%s.%s",d.getPdbId(),d.getRange());
			try {
				cache.getStructure(pdbRangeStr);
			} catch(Exception e) {
				logger.error("Can't get range "+pdbRangeStr,e);
				errors++;
				continue;
			}

			//All test passed
			logger.info("OK "+d.getDomainId());
		}
		return errors;
	}
}
