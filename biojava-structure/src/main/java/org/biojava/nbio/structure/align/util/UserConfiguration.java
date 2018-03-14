/*
 *                    PDB web development code
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
 *
 * Created on Jul 8, 2009
 * Created by ap3
 *
 */

package org.biojava.nbio.structure.align.util;

import org.biojava.nbio.structure.align.ce.StartupParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.biojava.nbio.core.util.XMLWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;


/** A container to persist config to the file system
 *
 * @author Andreas Prlic
 *
 */
public class UserConfiguration
{

	private static final Logger logger = LoggerFactory.getLogger(UserConfiguration.class);

	public static final String PDB_FORMAT   = "PDB";
	public static final String MMCIF_FORMAT = "mmCif";
	public static final String MMTF_FORMAT  = "mmtf";

	public static final String TMP_DIR = "java.io.tmpdir";

	public static final String PDB_DIR   = "PDB_DIR";
	public static final String PDB_CACHE_DIR = "PDB_CACHE_DIR";

	public static final String lineSplit = System.getProperty("file.separator");

	private String pdbFilePath;
	private String cacheFilePath;

	private FetchBehavior fetchBehavior;
	private ObsoleteBehavior obsoleteBehavior;

	private String fileFormat;

	private static AtomicBoolean warningShown = new AtomicBoolean(false);


	/**
	 * Default UserConfiguration:
	 * <ul>
	 * <li>split directory</li>
	 * <li>autofetch files</li>
	 * <li>default download location. This is the first specified of:
	 * 	<ol><li>{@value #PDB_DIR} system property (for instance, -D{@value #PDB_DIR}=/tmp)</li>
	 *   <li>{@value #PDB_DIR} environment variable</li>
	 *   <li>System temp directory (java.io.tmpdir property)</li>
	 *   </ol>
	 *   if the provided path is not a directory then
	 *   the system's temp directory is used. A non-writable path is allowed,
	 *   only a warning will be logged.
	 * </li>
	 * <li>default cache location. This is the first specified of:
	 * 	<ol><li>{@value #PDB_CACHE_DIR} system property (for instance, -D{@value #PDB_CACHE_DIR}=/tmp)</li>
	 *   <li>{@value #PDB_CACHE_DIR} environment variable</li>
	 *   <li>the value set for {@value #PDB_DIR}</li>
	 *   </ol>
	 *   if the provided path is not a directory or is not writable then
	 *   the system's temp directory is used.
	 * </li>
	 * </ul>
	 */
	public UserConfiguration(){
		fetchBehavior = FetchBehavior.DEFAULT;
		obsoleteBehavior = ObsoleteBehavior.DEFAULT;

		pdbFilePath = initPdbFilePath();
		// note that in initCacheFilePath, we set to the provided one (if readable) or to the same as pdbFilePath
		cacheFilePath = initCacheFilePath();

		fileFormat = MMTF_FORMAT;
	}

	private String initPdbFilePath() {

		String path = null;

		String propertyName = PDB_DIR;

		String userProvidedDir = System.getProperty(propertyName);

		if ( userProvidedDir != null && !userProvidedDir.trim().isEmpty()) {

			path = userProvidedDir;
			logger.debug("Read PDB dir from system property {}: {}", propertyName, path);
			File f = new File(path);
			if (!f.isDirectory()) {
				logger.warn(
						"Provided path {} (with system property {}) is not a directory. Using system's temp directory instead {}",
						path, propertyName, System.getProperty(TMP_DIR));
				path = System.getProperty(TMP_DIR);
			} else if (!f.canWrite()) {
				logger.warn(
						"Provided path {} (with system property {}) is not writable. Will not be able to write cached files.",
						path, propertyName);
				// we don't require the PDB_DIR to be writable, so that it can be used with a pre-rsynced dir
				// thus if not writable, we only warn and go ahead using it
			}


		} else {
			Map<String,String> env = System.getenv();

			if( env.containsKey(propertyName) && !env.get(propertyName).trim().isEmpty()) {
				path = env.get(propertyName);
				logger.debug("Read dir from environment variable {}: {}", propertyName, path);

				File f = new File(path);
				if (!f.isDirectory()) {
					logger.warn(
							"Provided path {} (with environment variable {}) is not a directory. Using system's temp directory instead {}",
							path, propertyName, System.getProperty(TMP_DIR));
					path = System.getProperty(TMP_DIR);
				} else if (!f.canWrite()) {
					logger.warn(
							"Provided path {} (with environment variable {}) is not writable. Will not be able to write cached files",
							path, propertyName);
					// we don't require the PDB_DIR to be writable, so that it can be used with a pre-rsynced dir
					// thus if not writable, we only warn and go ahead using it
				}

			} else {
				path = System.getProperty(TMP_DIR);

				if ( ! warningShown.get()) {

					logger.warn("Could not read dir from system property {} or environment variable {}, "
									+ "using system's temp directory {}",
							propertyName, propertyName, path);

					warningShown.set(true);
				}

				System.setProperty(propertyName,path);
			}
		}

		if ( ! path.endsWith(lineSplit) )
			path = path + lineSplit;

		return path;

	}

	private String initCacheFilePath() {

		String path = null;

		String propertyName = PDB_CACHE_DIR;

		String userProvidedDir = System.getProperty(propertyName);

		if ( userProvidedDir != null ) {

			path = userProvidedDir;
			logger.debug("Read cache dir from system property {}: {}", propertyName, path);
			File f = new File(path);
			if (!f.isDirectory()) {
				logger.warn(
						"Provided path {} (with system property {}) is not a directory. Using system's temp directory instead {}",
						path, propertyName, System.getProperty(TMP_DIR));
				path = System.getProperty(TMP_DIR);
			} else if (!f.canWrite()) {
				logger.warn(
						"Provided path {} (with system property {}) is not writable. Using system's temp directory instead {}",
						path, propertyName, System.getProperty(TMP_DIR));
				path = System.getProperty(TMP_DIR);
				System.setProperty(propertyName,path);
			}


		} else {
			Map<String,String> env = System.getenv();

			if( env.containsKey(propertyName)) {
				path = env.get(propertyName);
				logger.debug("Read dir from environment variable {}: {}", propertyName, path);

				File f = new File(path);
				if (!f.isDirectory()) {
					logger.warn(
							"Provided path {} (with environment variable {}) is not a directory. Using system's temp directory instead {}",
							path, propertyName, System.getProperty(TMP_DIR));
					path = System.getProperty(TMP_DIR);
				} else if (!f.canWrite()) {
					logger.warn(
							"Provided path {} (with environment variable {}) is not writable. Using system's temp directory instead {}",
							path, propertyName, System.getProperty(TMP_DIR));
					path = System.getProperty(TMP_DIR);
				}

			} else {
				// NOTE in case of not provided, then it is set to same as pdbFilePath
				// as PDB_DIR is not checked for being writable, we have to do that check here in case
				if (new File(pdbFilePath).canWrite()){
					path = pdbFilePath;
					logger.info("Could not read cache dir from system property {} or environment variable {}, "
							+ "using PDB directory instead {}",
							propertyName, propertyName, path);
					System.setProperty(propertyName,path);

				} else {
					path = System.getProperty(TMP_DIR);
					logger.warn("Could not read cache dir from system property {} or environment variable {}, "
							+ "and PDB directory {} is not writable. Using system's temp directory instead {}",
							propertyName, propertyName, pdbFilePath, path);
					System.setProperty(propertyName,path);

				}
			}
		}

		if ( ! path.endsWith(lineSplit) )
			path = path + lineSplit;

		return path;

	}

	public String getPdbFilePath()
	{
		return pdbFilePath;
	}

	public void setPdbFilePath(String pdbFilePath)
	{
		this.pdbFilePath = pdbFilePath;
	}

	public String getCacheFilePath() {
		return cacheFilePath;
	}

	public void setCacheFilePath(String cacheFilePath) {
		this.cacheFilePath = cacheFilePath;
	}

	public FetchBehavior getFetchBehavior() {
		return fetchBehavior;
	}

	public void setFetchBehavior(FetchBehavior fetchBehavior) {
		this.fetchBehavior = fetchBehavior;
	}

	public ObsoleteBehavior getObsoleteBehavior() {
		return obsoleteBehavior;
	}

	public void setObsoleteBehavior(ObsoleteBehavior obsoleteBehavior) {
		this.obsoleteBehavior = obsoleteBehavior;
	}

	/** convert Configuration to an XML file so it can be serialized
	 *
	 * @param pw
	 * @return XMLWriter
	 * @throws IOException
	 */
	public XMLWriter toXML(PrintWriter pw)
			throws IOException
			{

		XMLWriter     xw = new PrettyXMLWriter( pw);

		toXML(xw);
		return xw ;
			}


	/** convert Configuration to an XML file so it can be serialized
	 * add to an already existing xml file.
	 *
	 * @param xw the XML writer to use
	 * @return the writer again
	 * @throws IOException
	 * @see org.biojava.nbio.structure.align.webstart.ConfigXMLHandler
	 */

	public XMLWriter toXML(XMLWriter xw)
			throws IOException
			{
		xw.printRaw("<?xml version='1.0' standalone='no' ?>");
		//xw.printRaw("<!DOCTYPE " + XML_CONTENT_TYPE + " SYSTEM '" + XML_DTD + "' >");
		xw.openTag("JFatCatConfig");

		xw.openTag("PDBFILEPATH");
		// we don;t serialize the tempdir...
		String tempdir = System.getProperty(TMP_DIR);
		if (! pdbFilePath.equals(tempdir))
			xw.attribute("path", pdbFilePath);

		xw.attribute("fetchBehavior", fetchBehavior+"");
		xw.attribute("obsoleteBehavior", obsoleteBehavior+"");
		xw.attribute("fileFormat", fileFormat);
		xw.closeTag("PDBFILEPATH");

		xw.closeTag("JFatCatConfig");
		return xw ;

			}

	public static UserConfiguration fromStartupParams(StartupParameters params) {
		UserConfiguration config = new UserConfiguration();
		config.setPdbFilePath(params.getPdbFilePath());
		
		if(params.isAutoFetch()) {
			config.setFetchBehavior(FetchBehavior.DEFAULT);
		} else {
			config.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		}
		
		// TODO support MMCif Files
		config.setFileFormat(UserConfiguration.PDB_FORMAT);
		return config;
	}

	public void setFileFormat (String fileFormat){
		this.fileFormat = fileFormat;
	}

	public String getFileFormat()
	{
		return fileFormat;
	}






}
