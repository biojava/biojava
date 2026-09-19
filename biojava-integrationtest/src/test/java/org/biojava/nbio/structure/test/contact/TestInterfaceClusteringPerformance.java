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
package org.biojava.nbio.structure.test.contact;

import org.biojava.nbio.core.util.SingleLinkageClusterer;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.InterfaceFinder;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceCluster;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.StructureFiletype;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.assertEquals;

/**
 * Compares the performance of interface clustering with {@link StructureInterfaceList#clusterInterfaces}
 * (leader algorithm) against single linkage clustering with {@link SingleLinkageClusterer}, for an assembly with many interfaces.
 *
 * By default it is ignored, also by maven when the test is selected explicitly with -Dtest.
 * To execute it, run it from the IDE or remove the {@code @Ignore} annotation temporarily.
 */
public class TestInterfaceClusteringPerformance {

	/** An icosahedral capsid with 180 chains of one entity and ~700 interfaces. Larger cases like 5vf3 or 3j3q take very long */
	private static final String PDB_ID = "1gav";

	@Ignore("Performance test to be run manually")
	@Test
	public void testClusteringPerformance() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		cache.setFiletype(StructureFiletype.CIF);
		StructureIO.setAtomCache(cache);

		Structure assembly = StructureIO.getBiologicalAssembly(PDB_ID, 1, false);
		List<StructureInterface> list = new InterfaceFinder(assembly).getAllInterfaces().getList();
		System.out.printf("Found %d interfaces in assembly 1 of %s%n", list.size(), PDB_ID);

		double cutoff = StructureInterfaceList.DEFAULT_CONTACT_OVERLAP_SCORE_CLUSTER_CUTOFF;

		// first calculations of contact overlap scores initialise the contact sets (lazily): we do that before timing
		for (StructureInterface interf : list) {
			interf.getGroupContacts();
		}

		long start = System.currentTimeMillis();
		List<StructureInterfaceCluster> singletons = new ArrayList<>();
		for (StructureInterface interf : list) {
			StructureInterfaceCluster cluster = new StructureInterfaceCluster();
			cluster.addMember(interf);
			singletons.add(cluster);
		}
		List<StructureInterfaceCluster> leaderClusters = StructureInterfaceList.clusterInterfaces(singletons, StructureInterfaceList.ENTITY_ID_PAIR, cutoff);
		long end = System.currentTimeMillis();
		System.out.printf("%d clusters found. Time for leader clustering: %d ms%n", leaderClusters.size(), end - start);

		start = System.currentTimeMillis();
		double[][] matrix = new double[list.size()][list.size()];
		for (int i = 0; i < list.size(); i++) {
			for (int j = i + 1; j < list.size(); j++) {
				matrix[i][j] = Math.max(list.get(i).getContactOverlapScore(list.get(j), false), list.get(i).getContactOverlapScore(list.get(j), true));
			}
		}
		long endMatrix = System.currentTimeMillis();
		Map<Integer, Set<Integer>> singleLinkageClusters = new SingleLinkageClusterer(matrix, true).getClusters(cutoff);
		end = System.currentTimeMillis();
		System.out.printf("%d clusters found. Time for single linkage clustering: %d ms (all-vs-all scores %d ms, clustering %d ms)%n",
				singleLinkageClusters.size(), end - start, endMatrix - start, end - endMatrix);

		assertEquals(list.size(), leaderClusters.stream().mapToInt(c -> c.getMembers().size()).sum());
		assertEquals(list.size(), singleLinkageClusters.values().stream().mapToInt(Set::size).sum());
	}
}
