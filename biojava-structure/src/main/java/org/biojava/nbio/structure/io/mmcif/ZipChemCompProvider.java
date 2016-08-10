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
package org.biojava.nbio.structure.io.mmcif;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** This chemical component provider retrieves and caches chemical component definition files from a
 * zip archive specified in its construction.  If the archive does not contain the record, an attempt is
 * made to download it using DownloadChemCompProvider. The downloaded file is then added to the archive.
 *
 * The class is thread-safe and the same ZipChemCompProvider should be used by all threads to prevent
 * simultaneous read or write to the zip archive.  A zip archive will be created if missing.
 *
 * @author edlunde
 * @author larsonm
 * @since 12/05/12
 * updated 3/5/2016 for Java 7 ZipFileSystem
 */
public class ZipChemCompProvider implements ChemCompProvider{
	private static final Logger s_logger = LoggerFactory.getLogger(ZipChemCompProvider.class);

	private final Path m_tempDir;  // Base path where $m_zipRootDir/ will be downloaded to.
	private final Path m_zipRootDir;
	private final Path m_zipFile;
	private final DownloadChemCompProvider m_dlProvider;

	private boolean m_removeCif;

	// Missing IDs from library that cannot be download added here to prevent delays.
	private Set<String> unavailable = new HashSet<String>();

	/**
	 * ZipChemCompProvider is a Chemical Component provider that stores chemical components
	 * in a zip archive.  Missing chemical components are downloaded and appended to the
	 * archive.  If non-existent a new zip archive will be created.
	 *
	 * @param chemicalComponentDictionaryFile : path to zip archive for chemical components.
	 * @param tempDir : path for temporary directory, (null) defaults to path in property "java.io.tmpdir".
	 * @throws IOException
	 */
	public ZipChemCompProvider(String chemicalComponentDictionaryFile, String tempDir) throws IOException {
		this.m_zipFile = Paths.get(chemicalComponentDictionaryFile);

		// Use a default temporary directory if not passed a value.
		if (tempDir == null || tempDir.equals("")) {
			this.m_tempDir = Paths.get(System.getProperty("java.io.tmpdir"));
		} else {
			this.m_tempDir = Paths.get(tempDir);
		}

		this.m_zipRootDir = Paths.get("chemcomp");

		// Setup an instance of the download chemcomp provider.
		this.m_dlProvider = new DownloadChemCompProvider(m_tempDir.toString());
		this.m_removeCif = true;
		initializeZip();
	}

	// See comments in addToZipFileSystem for why initialization is required with
	// ZipFileSystems - due to URI issues in Java7.
	private void initializeZip() throws IOException {
		s_logger.info("Using chemical component dictionary: " + m_zipFile.toString());
		final File f = m_zipFile.toFile();
		if (!f.exists()) {
			s_logger.info("Creating missing zip archive: " + m_zipFile.toString());
			FileOutputStream fo = new FileOutputStream(f);
			ZipOutputStream zip = new ZipOutputStream(new BufferedOutputStream(fo));
			try {
				zip.putNextEntry(new ZipEntry("chemcomp/"));
				zip.closeEntry();
			} finally {
				zip.close();
			}
		}
	}

	/**
	 * Remove downloaded .cif.gz after adding to zip archive?
	 * Default is true.
	 * @param doRemove
	 */
	public void setRemoveCif(boolean doRemove) {
		m_removeCif = doRemove;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.io.mmcif.ChemCompProvider#getChemComp(java.lang.String)
	 *
	 * @param recordName : three letter PDB name for a residue
	 * @return ChemComp from .zip or ChemComp from repository.  Will return empty ChemComp when unable to find a residue and will return null if not provided a valid recordName.
	 */
	@Override
	public ChemComp getChemComp(String recordName) {
		if (null == recordName) return null;

		// handle non-existent ChemComp codes and do not repeatedly attempt to add these.
		for (String str : unavailable) {
			if (recordName.equals(str)) return getEmptyChemComp(recordName);
		}

		// Try to pull from zip, if fail then download.
		ChemComp cc = getFromZip(recordName);
		if (cc == null) {
			s_logger.info("File "+recordName+" not found in archive. Attempting download from PDB.");
			cc = downloadAndAdd(recordName);
		}

		// If a null record or an empty chemcomp, return a default ChemComp and blacklist.
		if (cc == null || (null == cc.getName() && cc.getAtoms().size() == 0)) {
			s_logger.info("Unable to find or download " + recordName + " - excluding from future searches.");
			unavailable.add(recordName);
			return getEmptyChemComp(recordName);
		}
		return cc;
	}

	/** Use DownloadChemCompProvider to grab a gzipped cif record from the PDB.
	 *  Zip all downloaded cif.gz files into the dictionary.
	 *
	 * @param recordName is the three-letter chemical component code (i.e. residue name).
	 * @return ChemComp matching recordName
	 */
	private ChemComp downloadAndAdd(String recordName){
		final ChemComp cc = m_dlProvider.getChemComp(recordName);

		// final File [] files = finder(m_tempDir.resolve("chemcomp").toString(), "cif.gz");
		final File [] files = new File[1];
		Path cif = m_tempDir.resolve("chemcomp").resolve(recordName + ".cif.gz");
		files[0] = cif.toFile();
		if (files[0] != null) {
			addToZipFileSystem(m_zipFile, files, m_zipRootDir);
			if (m_removeCif) for (File f : files) f.delete();
		}
		return cc;
	}

	/**
	 * Cleanup chemical component (.cif.gz) files downloaded to tmpdir.
	 * @param tempdir : path to temporary directory for chemical components
	 */
	public static void purgeTempFiles(String tempdir) {
		if (tempdir == null) return;

		s_logger.info("Removing: "+tempdir);
		Path dlPath = Paths.get(tempdir).resolve("chemcomp");
		File[] chemCompOutFiles = finder(dlPath.toString(), "cif.gz");
		if (null != chemCompOutFiles) for (File f : chemCompOutFiles) f.delete();
		dlPath.toFile().delete();
	}

	/**
	 * Return an empty ChemComp group for a three-letter resName.
	 * @param resName
	 * @return
	 */
	private ChemComp getEmptyChemComp(String resName){
		String pdbName = ""; // Empty string is default
		if (null != resName && resName.length() >= 3) {
			pdbName = resName.substring(0,3);
		}
		final ChemComp comp = new ChemComp();
		comp.setOne_letter_code("?");
		comp.setThree_letter_code(pdbName);
		comp.setPolymerType(PolymerType.unknown);
		comp.setResidueType(ResidueType.atomn);
		return comp;
	}

	/**
	 * Return File(s) in dirName that match suffix.
	 * @param dirName
	 * @param suffix
	 * @return
	 */
	static private File[] finder( String dirName, final String suffix){
		if (null == dirName || null == suffix) {
			return null;
		}

		final File dir = new File(dirName);
		return dir.listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String filename)
			{ return filename.endsWith(suffix); }
		} );
	}

	/**
	 * This is synchronized, along with addToFileSystem to prevent simulatenous reading/writing.
	 * @param recordName to find in zipfile.
	 * @return ChemComp if found or null if missing.
	 */
	private synchronized ChemComp getFromZip(String recordName) {
		ChemComp cc = null;
		if (!m_zipFile.toFile().exists()) return cc;
		final String filename = "chemcomp/" + recordName+".cif.gz";

		// try with resources block to read from the filesystem.
		try (FileSystem fs = FileSystems.newFileSystem(m_zipFile, null)) {
			Path cif = fs.getPath(filename);

			if (Files.exists(cif)) {
				final InputStream zipStream = Files.newInputStream(cif);
				final InputStream inputStream = new GZIPInputStream(zipStream);
				s_logger.debug("reading " + recordName + " from " + m_zipFile);
				final MMcifParser parser = new SimpleMMcifParser();
				final ChemCompConsumer consumer = new ChemCompConsumer();
				parser.addMMcifConsumer(consumer);
				parser.parse(inputStream);
				inputStream.close();

				final ChemicalComponentDictionary dict = consumer.getDictionary();
				cc = dict.getChemComp(recordName);
			}
		} catch (IOException e) {
			s_logger.error("Unable to read from zip file : " + e.getMessage());
		}

		return cc;
	}

	/**
	 * Add an array of files to a zip archive.
	 * Synchronized to prevent simultaneous reading/writing.
	 *
	 * @param zipFile is a destination zip archive
	 * @param files is an array of files to be added
	 * @param pathWithinArchive is the path within the archive to add files to
	 * @return true if successfully appended these files.
	 */
	private synchronized boolean addToZipFileSystem(Path zipFile, File[] files, Path pathWithinArchive) {
		boolean ret = false;

		/* URIs in Java 7 cannot have spaces, must use Path instead
		 * and so, cannot use the properties map to describe need to create
		 * a new zip archive.  ZipChemCompProvider.initilizeZip to creates the
		 * missing zip file */

		/*
		// convert the filename to a URI
		String uriString = "jar:file:" + zipFile.toUri().getPath();
		final URI uri = URI.create(uriString);

		// if filesystem doesn't exist, create one.
		final Map<String, String> env = new HashMap<>();
		// Create a new zip if one isn't present.
		if (!zipFile.toFile().exists()) {
			System.out.println("Need to create " + zipFile.toString());
		}
		env.put("create", String.valueOf(!zipFile.toFile().exists()));
		// Specify the encoding as UTF -8
		env.put("encoding", "UTF-8");
		*/

		// Copy in each file.
		try (FileSystem zipfs = FileSystems.newFileSystem(zipFile, null)) {
			Files.createDirectories(pathWithinArchive);
			for (File f : files) {
				if (!f.isDirectory() && f.exists()) {
					Path externalFile = f.toPath();
					Path pathInZipFile = zipfs.getPath(pathWithinArchive.resolve(f.getName()).toString());
					Files.copy(externalFile, pathInZipFile,
							StandardCopyOption.REPLACE_EXISTING);
				}
			}
			ret = true;
		} catch (IOException ex) {
			s_logger.error("Unable to add entries to Chemical Component zip archive : " + ex.getMessage());
			ret = false;
		}
		return ret;
	}
}
