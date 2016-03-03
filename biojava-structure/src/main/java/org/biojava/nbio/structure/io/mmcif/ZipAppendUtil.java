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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URI;
import java.nio.file.*;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import java.util.zip.*;

/**
 * Singleton helper class to add entries to zip archives.
 * @since Feb 8, 2013
 * @author edlunde
 * @author larsonm
 *
 */
public class ZipAppendUtil {
	private static final Logger s_logger = Logger.getLogger( ZipAppendUtil.class.getPackage().getName() );
	private static ZipAppendUtil s_instance = null; // Singleton instance
	private boolean m_busy = false;

	/**
	 * Get a singleton instance.
	 * @return
	 */
	public static ZipAppendUtil getInstance() {
		if (s_instance == null) {
			s_instance = new ZipAppendUtil();
		}
		return s_instance;
	}

	/**
	 * Add an array of files to a zip archive.
	 *
	 * @param zipFile is a destination zip archive
	 * @param files is an array of files to be added
	 * @param pathWithinArchive is the path within the archive to add files to
     * @return true if successfully appended these files.
     */
	public boolean addToZipFileSystem(Path zipFile, File[] files, Path pathWithinArchive) {
		boolean ret = false;
		if (m_busy) return ret; // If another thread is adding files, return immediately.
		m_busy = true;

		// Copy in each file.
		try (FileSystem zipfs = createZipFileSystem(zipFile.toString())) {
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
			s_logger.severe("Unable to add entries to Chemical Component zip archive : " + ex.getMessage());
			ret = false;
		}
		m_busy = false;
		return ret;
	}

	/**
	 * Returns a zip file system
	 * @param zipFilename to construct the file system from
	 * @return a zip file system
	 * @throws IOException
	 */
	private static FileSystem createZipFileSystem(String zipFilename) throws IOException {
		// convert the filename to a URI
		final Path path = Paths.get(zipFilename);
		final URI uri = URI.create("jar:file:" + path.toUri().getPath());

		final Map<String, String> env = new HashMap<>();
		return FileSystems.newFileSystem(uri, env);
	}

	/**
	 * Create a copy for a file.
	 * May be necessary to use a copy of the zip archive when appending.
	 *
	 * @param src
	 * @param dst
	 * @throws IOException
	 */
	private void copyFile(File src, File dst) throws IOException {
	    InputStream is = null;
	    OutputStream os = null;
	    try {
	        is = new FileInputStream(src);
	        os = new FileOutputStream(dst);
	        byte[] buffer = new byte[1024];
	        int length;
	        while ((length = is.read(buffer)) > 0) {
	            os.write(buffer, 0, length);
	        }
	    } finally {
	        if (is != null) is.close();
	        if (os != null) os.close();
	    }
	}
	
	/* Prevent direct access to the constructor*/
	private ZipAppendUtil() {
		super();
	}

	/**
	 * @return true if zipFile is in use
     */
	public boolean isBusy() {
		return m_busy;
	}
}
