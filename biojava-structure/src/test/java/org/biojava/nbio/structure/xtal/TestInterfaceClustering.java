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
package org.biojava.nbio.structure.xtal;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.io.StructureFiletype;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.asa.GroupAsa;
import org.biojava.nbio.structure.contact.InterfaceFinder;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.contact.StructureInterface;
import org.biojava.nbio.structure.contact.StructureInterfaceCluster;
import org.biojava.nbio.structure.contact.StructureInterfaceList;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.junit.Test;

import javax.vecmath.Matrix4d;
import java.util.function.Function;

public class TestInterfaceClustering {

	@Test
	public void test3DDO() throws IOException, StructureException {

		// 3DDO is special in that it contains 6 chains in 1 entity, all of them with different residue numbering

		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		cache.setFiletype(StructureFiletype.CIF);

		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("3DDO");

		CrystalBuilder cb = new CrystalBuilder(s);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		interfaces.calcAsas(100, 1, 0);
		interfaces.removeInterfacesBelowArea();

		List<StructureInterfaceCluster> clusters = interfaces.getClusters();

		// 22 if below 35A2 interfaces are filtered
		assertEquals(22,interfaces.size());

		// we simply want to test that some interfaces cluster together, for this entry
		// it is problematic because of different residue numbering between different chains of same entity
		assertTrue("Expected fewer than 22 interfaces (some interfaces should cluster together)",clusters.size()<22);

		// first 2 clusters are of size 3
		assertEquals("Cluster 1 should have 3 members",3,clusters.get(0).getMembers().size());
		assertEquals("Cluster 2 should have 3 members",3,clusters.get(1).getMembers().size());

		// detection of isologous test: first 3 interfaces should be isologous

		assertTrue("Interface 1 should be isologous",interfaces.get(1).isIsologous());
		assertTrue("Interface 2 should be isologous",interfaces.get(2).isIsologous());
		assertTrue("Interface 3 should be isologous",interfaces.get(3).isIsologous());



	}

	/**
	 * Test for NCS clustering in viral capsid structures that contain NCS operators.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test1AUY() throws IOException, StructureException {

		// 1AUY is a viral capsid with NCS ops

		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		cache.setFiletype(StructureFiletype.CIF);

		StructureIO.setAtomCache(cache);

		// 3vbr would be an example of capsids with several chains
		Structure s = StructureIO.getStructure("1auy");

		Map<String,String> chainOrigNames = new HashMap<>();
		Map<String,Matrix4d> chainNcsOps = new HashMap<>();
		CrystalBuilder.expandNcsOps(s,chainOrigNames,chainNcsOps);
		CrystalBuilder cb = new CrystalBuilder(s, chainOrigNames, chainNcsOps);

		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);

		List<StructureInterfaceCluster> clusters = interfaces.getClusters();

		assertNotNull(clusters);

		assertTrue(clusters.size()<=interfaces.size());

		interfaces.calcAsas(100, 1, 0);

		// after calculating ASAs we should have ids for all interfaces
		for (StructureInterface interf : interfaces) {
			assertTrue(interf.getId()>0);
		}


		int numInterfacesShouldbeKept = 0;

		List<StructureInterfaceCluster> ncsClusterShouldbeKept = new ArrayList<>();
		for (StructureInterfaceCluster ncsCluster : interfaces.getClustersNcs()) {
			if (ncsCluster.getMembers().get(0).getTotalArea()>=StructureInterfaceList.DEFAULT_MINIMUM_INTERFACE_AREA) {
				//System.out.println("NCS cluster is above cutoff area and has "+ncsCluster.getMembers().size()+ " members");
				ncsClusterShouldbeKept.add(ncsCluster);
				numInterfacesShouldbeKept += ncsCluster.getMembers().size();
			}
		}

		clusters = interfaces.getClusters();

		assertNotNull(clusters);

		assertTrue(clusters.size()<=interfaces.size());

		interfaces.removeInterfacesBelowArea();

		assertNotNull(interfaces.getClustersNcs());

		// making sure that removeInterfacesBelowArea does not throw away the members for which area wasn't calculated
		for (StructureInterfaceCluster ncsCluster : ncsClusterShouldbeKept) {
			assertTrue(interfaces.getClustersNcs().contains(ncsCluster));
		}

		assertEquals(numInterfacesShouldbeKept, interfaces.size());

		clusters = interfaces.getClusters();

		assertNotNull(clusters);

		assertTrue(clusters.size()<=interfaces.size());

		for (StructureInterface interf : interfaces) {
			GroupAsa groupAsa = interf.getFirstGroupAsas().values().iterator().next();
			String expected = interf.getMoleculeIds().getFirst();
			String actual = groupAsa.getGroup().getChain().getName();
			// in 1auy this works always since there's only 1 chain. But it is useful in testing cases like 3vbr with serveral chains
			assertEquals(expected.charAt(0), actual.charAt(0));
		}
	}


	@Test
	public void test3C5FWithSeqresPdb() throws IOException, StructureException {

		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/3c5f_raw.pdb.gz"));
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(s);

		assertEquals(8, s.getPolyChains().size());

		// 1 protein, 3 nucleotide chains, 1 NA nonpoly chain, 1 water: 6 entities
		assertEquals(6, s.getEntityInfos().size());

		CrystalBuilder cb = new CrystalBuilder(s);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		interfaces.calcAsas(100, 1, 0);
		interfaces.removeInterfacesBelowArea();

		List<StructureInterfaceCluster> clusters = interfaces.getClusters();

		// 23 if below 35A2 interfaces are filtered
		assertEquals(23,interfaces.size());

		// we simply want to test that some interfaces cluster together
		assertTrue("Expected fewer than 23 interfaces (some interfaces should cluster together)",clusters.size()<23);

		// third cluster (index 2) is of size 2
		assertEquals("Cluster 3 should have 2 members",2,clusters.get(2).getMembers().size());

		assertTrue("Interface 3 should be isologous",interfaces.get(3).isIsologous());


	}

	// This doesn't work yet, since for raw files without a SEQRES, the seqres groups are not populated. Instead
	// in that case Compound.getAlignedResIndex() returns residue numbers as given (without insertion codes) and
	// thus in general residues will not be correctly aligned between different chains of same entity. This breaks
	// cases like 3ddo (with no SEQRES records) where residue numbering is different in every chain of the one entity.
	// Then contact overlap calculation will be wrong and interface clustering won't work.
	// see https://github.com/eppic-team/eppic/issues/39
	// See also TestCompoundResIndexMapping
	//@Test
	public void test3DDONoSeqresPdb() throws IOException, StructureException {

		// 3ddo contains 6 chains in 1 entity, with residue numbering completely different in each of the chains

		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/3ddo_raw_noseqres.pdb.gz"));
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(s);

		assertEquals(6, s.getChains().size());

		assertEquals(1, s.getEntityInfos().size());

		CrystalBuilder cb = new CrystalBuilder(s);
		StructureInterfaceList interfaces = cb.getUniqueInterfaces(5.5);
		interfaces.calcAsas(100, 1, 0);
		interfaces.removeInterfacesBelowArea();

		List<StructureInterfaceCluster> clusters = interfaces.getClusters();

		// 22 if below 35A2 interfaces are filtered
		assertEquals(22,interfaces.size());

		// we simply want to test that some interfaces cluster together, for this entry
		// it is problematic because of different residue numbering between different chains of same entity
		assertTrue("Expected fewer than 22 interfaces (some interfaces should cluster together)",clusters.size()<22);

		// first 2 clusters are of size 3
		assertEquals("Cluster 1 should have 3 members",3,clusters.get(0).getMembers().size());
		assertEquals("Cluster 2 should have 3 members",3,clusters.get(1).getMembers().size());

		// detection of isologous test: first 3 interfaces should be isologous

		assertTrue("Interface 1 should be isologous",interfaces.get(1).isIsologous());
		assertTrue("Interface 2 should be isologous",interfaces.get(2).isIsologous());
		assertTrue("Interface 3 should be isologous",interfaces.get(3).isIsologous());



	}

	/**
	 * Two-step clustering of the interfaces of an assembly: first geometrically identical interfaces (same asym ids),
	 * then similar interfaces between same entities, starting from the clusters of the first step.
	 */
	@Test
	public void testClusterInterfacesTwoSteps() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		cache.setFiletype(StructureFiletype.CIF);
		StructureIO.setAtomCache(cache);

		// 6cco assembly 1: D3 with 2 entities, 12 interfaces
		Structure assembly = StructureIO.getBiologicalAssembly("6cco", 1, false);
		List<StructureInterface> interfaces = new InterfaceFinder(assembly).getAllInterfaces().getList();
		assertEquals(12, interfaces.size());

		Function<StructureInterface, Pair<String>> asymIdPair = interf -> new Pair<>(
				stripOperators(interf.getParentChains().getFirst().getId()),
				stripOperators(interf.getParentChains().getSecond().getId()));

		List<StructureInterfaceCluster> singletons = new ArrayList<>();
		for (StructureInterface interf : interfaces) {
			StructureInterfaceCluster cluster = new StructureInterfaceCluster();
			cluster.addMember(interf);
			singletons.add(cluster);
		}

		List<StructureInterfaceCluster> identicalClusters = StructureInterfaceList.clusterInterfaces(singletons, asymIdPair, 0.99);
		assertEquals(4, identicalClusters.size());
		for (StructureInterfaceCluster cluster : identicalClusters) {
			assertEquals(3, cluster.getMembers().size());
			cluster.getMembers().forEach(m -> assertSame(cluster, m.getCluster()));
		}
		// input is not modified
		singletons.forEach(c -> assertEquals(1, c.getMembers().size()));

		List<StructureInterfaceCluster> entityClusters = StructureInterfaceList.clusterInterfaces(identicalClusters, StructureInterfaceList.ENTITY_ID_PAIR, StructureInterfaceList.DEFAULT_CONTACT_OVERLAP_SCORE_CLUSTER_CUTOFF);
		assertEquals(3, entityClusters.size());
		int total = 0;
		for (StructureInterfaceCluster cluster : entityClusters) {
			// representatives are preserved from the first step
			assertTrue(identicalClusters.stream().anyMatch(c -> c.getMembers().get(0) == cluster.getMembers().get(0)));
			cluster.getMembers().forEach(m -> assertSame(cluster, m.getCluster()));
			total += cluster.getMembers().size();
		}
		assertEquals(12, total);
		identicalClusters.forEach(c -> assertEquals(3, c.getMembers().size()));
	}

	private static String stripOperators(String chainId) {
		return chainId.substring(0, chainId.indexOf(BiologicalAssemblyBuilder.SYM_CHAIN_ID_SEPARATOR));
	}

}
