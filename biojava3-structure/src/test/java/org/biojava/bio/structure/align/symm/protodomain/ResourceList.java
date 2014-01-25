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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure.align.symm.protodomain;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.symm.CeSymm;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.custommonkey.xmlunit.DetailedDiff;
import org.custommonkey.xmlunit.Diff;
import org.custommonkey.xmlunit.Difference;
import org.custommonkey.xmlunit.DifferenceListener;
import org.custommonkey.xmlunit.XMLUnit;
import org.custommonkey.xmlunit.examples.RecursiveElementNameAndTextQualifier;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;

/**
 * A singleton handler for testing resources, particularly alignment results but also other files. For every file, it's nice way to avoid hardcoding the resource path. Loads an alignment result if
 * they're available; otherwise, performs the alignment, saves it, and then loads it. When using ResourceList, avoid using hardcoded paths, and use {@link #getCache()} when an AtomCache is needed.
 * Uses {@link BerkeleyScopInstallation}, and as a <em>side effect calls {@link ScopFactory#setScopDatabase(org.biojava.bio.structure.scop.ScopDatabase)}</em>.
 * 
 * @author dmyerstu
 */
public class ResourceList {

	public static class ElementTextIgnoringDifferenceListener implements DifferenceListener {

		private String[] ignoredNames;

		public ElementTextIgnoringDifferenceListener(String... ignoredNames) {
			this.ignoredNames = ignoredNames;
		}

		public ElementTextIgnoringDifferenceListener(List<String> ignoredNames) {
			this.ignoredNames = new String[ignoredNames.size()];
			for (int i = 0; i < ignoredNames.size(); i++) {
				this.ignoredNames[i] = ignoredNames.get(i);
			}
		}

		@Override
		public int differenceFound(Difference difference) {
			Node controlNode = difference.getControlNodeDetail().getNode();
			String name = null;
			if (controlNode != null) {
				name = controlNode.getParentNode().getNodeName();
			}
			for (String ignoredName : ignoredNames) {
				if (ignoredName.equals(name)) {
					return RETURN_IGNORE_DIFFERENCE_NODES_IDENTICAL;
				}
			}
			return RETURN_ACCEPT_DIFFERENCE;
		}

		@Override
		public void skippedComparison(Node control, Node test) {
		}

	}

	/**
	 * Provides the main resource directory path, and full paths of alignment result files.
	 * 
	 * @author dmyerstu
	 * 
	 */
	public static abstract class NameProvider {

		/**
		 * @return A simple NameProvider using {@link NameProvider#DEFAULT_PDB_DIR} and the format {@code algorithm/nameA,nameB.xml} for alignments.
		 */
		public static NameProvider defaultNameProvider() {
			return new NameProvider() {
				@Override
				public String getNameForAlignment(String algorithmName, String nameA, String nameB) {
					return DEFAULT_DIR + algorithmName + "/" + nameA + "," + nameB + ".xml";
				}

				@Override
				public String getResourceDir() {
					return DEFAULT_DIR;
				}
			};
		}

		/**
		 * @return True iff the alignment by {@code algorithmName} of {@code nameA} against {@code nameB}, <em>in that order</em>. Not generally needed, since ResourceList tries to make this
		 *         transparent.
		 */
		public boolean alignmentExists(String algorithmName, String nameA, String nameB) {
			return getAlignmentFile(algorithmName, nameA, nameB).exists();
		}

		/**
		 * @return The File containing the XML-encoded AFPChain corresponding to the alignment by {@code algorithmName} of {@code nameA} against {@code nameB}, <em>in that order</em>.
		 */
		public File getAlignmentFile(String algorithmName, String nameA, String nameB) {
			File file = new File(getNameForAlignment(algorithmName, nameA, nameB));
			file.getParentFile().mkdirs();
			return file;
		}

		/**
		 * @return The file path containing the XML-encoded AFPChain corresponding to the alignment by {@code algorithmName} of {@code nameA} against {@code nameB}, <em>in that order</em>.
		 */
		public abstract String getNameForAlignment(String algorithmName, String nameA, String nameB);

		/**
		 * @return The main directory containing testing resources.
		 */
		public abstract String getResourceDir();
	}

	class AlignmentPair {
		private final String name1;
		private final String name2;

		public AlignmentPair(String name1, String name2) {
			this.name1 = name1;
			this.name2 = name2;
		}

		public String getName1() {
			return name1;
		}

		public String getName2() {
			return name2;
		}
	}

	public static final String DEFAULT_DIR = "src/test/resources/";

	/**
	 * Null means use the AtomCache default.
	 */
	public static final String DEFAULT_PDB_DIR = null;

	private static ResourceList singleton;

	static {
		StructureAlignmentFactory.addAlgorithm(new CeMain());
		StructureAlignmentFactory.addAlgorithm(new CeCPMain());
	}

	// notice the side effects here
	static {
		XMLUnit.setIgnoreWhitespace(true);
		XMLUnit.setIgnoreComments(true);
		XMLUnit.setIgnoreAttributeOrder(true);
	}

	/**
	 * Returns true if and only if the two XML files are <em>similar</em>; that is, the contain the same elements and attributes regardless of order.
	 */
	public static boolean compareXml(File expectedFile, File actualFile, DifferenceListener listener) {
		try {
			FileReader expectedFr = new FileReader(expectedFile);
			FileReader actualFr = new FileReader(actualFile);
			Diff diff = new Diff(expectedFr, actualFr);
			if (listener != null) diff.overrideDifferenceListener(listener);
			// ignore order
			// look at element, id, and weight (weight is a nested element)
			diff.overrideElementQualifier(new RecursiveElementNameAndTextQualifier());
			final boolean isSimilar = diff.similar();
			if (!isSimilar)
				printDetailedDiff(diff, System.err);
			expectedFr.close();
			actualFr.close();
			return isSimilar;
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (SAXException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * @return The singleton ResourceList.
	 */
	public static ResourceList get() {
		return singleton;
	}

	public static void printDetailedDiff(Diff diff, PrintStream ps) {
		DetailedDiff detDiff = new DetailedDiff(diff);
		for (Object object : detDiff.getAllDifferences()) {
			Difference difference = (Difference) object;
			ps.println(difference);
		}
	}

	/**
	 * Sets the singleton ResourceList to a new ResourceList with the given NameProvider and PDB directory. As a
	 * <em>side effect calls {@link ScopFactory#setScopDatabase(org.biojava.bio.structure.scop.ScopDatabase)}</em>.
	 */
	public static void set(NameProvider nameProvider, String pdbDir) {
		ResourceList.singleton = new ResourceList(nameProvider, pdbDir);
		ScopDatabase scop = ScopFactory.getSCOP();
		if (!scop.getClass().getName().equals(BerkeleyScopInstallation.class.getName())) { // for efficiency
			ScopFactory.setScopDatabase(new BerkeleyScopInstallation()); // ScopDatabase is too hard to mock well
		}
	}

	private static AFPChain fromXML(File file, String nameA, String nameB, Atom[] ca1, Atom[] ca2) throws IOException,
			StructureException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		StringBuilder sb = new StringBuilder();
		String line = "";
		while ((line = br.readLine()) != null)
			sb.append(line);
		br.close();
		return AFPChainXMLParser.fromXML(sb.toString(), ca1, ca2);
	}

	private AtomCache cache;

	private final ReadWriteLock lock = new ReentrantReadWriteLock();

	private NameProvider nameProvider;

	private ResourceList(NameProvider nameProvider, String pdbDir) {
		this.nameProvider = nameProvider;
		if (pdbDir == null) {
			cache = new AtomCache();
		} else {
			cache = new AtomCache(pdbDir, false);
		}
	}

	public Atom[] getAtoms(String name) {
		try {
			return cache.getAtoms(name);
		} catch (IOException e) {
			throw new ResourceException(e);
		} catch (StructureException e) {
			throw new ResourceException(e);
		}
	}

	/**
	 * Prefer using this to creating a new AtomCache, to avoid concurrency issues
	 * 
	 * @return
	 */
	public AtomCache getCache() {
		return cache;
	}

	public NameProvider getNameProvider() {
		return nameProvider;
	}

	public Structure getStructure(String name) {
		try {
			return cache.getStructure(name);
		} catch (IOException e) {
			throw new ResourceException(e);
		} catch (StructureException e) {
			throw new ResourceException(e);
		}
	}

	/**
	 * Loads the alignment requested, performing it iff necessary.
	 */
	public AFPChainAndAtoms load(File file, String nameA, String nameB) {
		System.out.println("Loading: " + file.getPath());
		try {
			Atom[] ca1 = cache.getAtoms(nameA);
			Atom[] ca2 = cache.getAtoms(nameB);
			AFPChain afpChain = fromXML(file, nameA, nameB, ca1, ca2);
			return new AFPChainAndAtoms(afpChain, ca1, ca2);
		} catch (Exception e) {
			throw new ResourceException("Could load the alignment.", e);
		}
	}

	/**
	 * Loads the alignment requested, performing it iff necessary.
	 */
	public AFPChainAndAtoms load(String algName, AlignmentPair names) {
		try {
			return load(StructureAlignmentFactory.getAlgorithm(algName), names);
		} catch (Exception e) {
			throw new ResourceException("Did not find the algorithm.", e);
		}
	}

	/**
	 * Loads the alignment requested, performing it iff necessary.
	 */
	public AFPChainAndAtoms load(StructureAlignment alg, AlignmentPair names) {
		AFPChainAndAtoms acaa;
		if (!nameProvider.alignmentExists(alg.getAlgorithmName(), names.getName1(), names.getName2())) {
			lock.writeLock().lock();
			acaa = put(alg, names.getName1(), names.getName2());
			lock.writeLock().unlock();
		} else {
			lock.readLock().lock();
			acaa = load(nameProvider.getAlignmentFile(alg.getAlgorithmName(), names.getName1(), names.getName2()),
					names.getName1(), names.getName2());
			lock.readLock().unlock();
		}
		return acaa;
	}

	public AFPChainAndAtoms load(StructureAlignment alg, String nameA, String nameB) {
		return load(alg, new AlignmentPair(nameA, nameB));
	}

	/**
	 * Loads the alignment with {@link CeMain} requested, performing it iff necessary.
	 */
	public AFPChainAndAtoms loadSim(String name1, String name2) {
		return load(new CeMain(), name1, name2);
	}

	/**
	 * Loads the alignment with {@link CeSymm} requested, performing it iff necessary.
	 */
	public AFPChainAndAtoms loadSymm(String name) {
		return load(new CeSymm(), name, name);
	}

	public File openFile(String filename) {
		return new File(nameProvider.getResourceDir() + filename);
	}

	public String openFileAsString(File file) {
		try {
			StringBuilder sb = new StringBuilder();
			BufferedReader br = openReader(file.getPath());
			String line = "";
			while ((line = br.readLine()) != null) {
				sb.append(line + "\n");
			}
			br.close();
			return sb.toString();
		} catch (IOException e) {
			throw new ResourceException(e);
		}
	}

	public String openFileAsString(String filename) {
		return openFileAsString(new File(filename));
	}

	public BufferedReader openReader(File file) {
		return openReader(file.getPath());
	}

	/**
	 * Opens the file as a {@link BufferedReader}.
	 */
	public BufferedReader openReader(String filename) {
		File file = openFile(filename);
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e) {
			throw new ResourceException(e);
		}
		return br;
	}

	/**
	 * Opens the file as a {@link FileInputStream}.
	 */
	public FileInputStream openStream(String filename) {
		File file = openFile(filename);
		FileInputStream fis;
		try {
			fis = new FileInputStream(file);
		} catch (FileNotFoundException e) {
			throw new ResourceException(e);
		}
		return fis;
	}

	public void setCache(AtomCache cache) {
		this.cache = cache;
	}

	private AFPChainAndAtoms put(StructureAlignment alg, String nameA, String nameB) {
		System.out.println("Putting: " + nameProvider.getAlignmentFile(alg.getAlgorithmName(), nameA, nameB).getPath());
		try {
			Atom[] ca1 = cache.getAtoms(nameA);
			Atom[] ca2 = cache.getAtoms(nameB);
			AFPChain afpChain = alg.align(ca1, ca2);
			afpChain.setName1(nameA);
			afpChain.setName2(nameB);
			String s = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n";
			s += AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
			BufferedWriter br = new BufferedWriter(new FileWriter(nameProvider.getAlignmentFile(alg.getAlgorithmName(),
					nameA, nameB)));
			br.write(s);
			br.close();
			afpChain.setName1(nameA);
			afpChain.setName2(nameB);
			return new AFPChainAndAtoms(afpChain, ca1, ca2);
		} catch (Exception e) {
			throw new ResourceException("Could save the alignment result.", e);
		}
	}

}
