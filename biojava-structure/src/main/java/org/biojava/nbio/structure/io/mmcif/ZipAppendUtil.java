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
import java.util.logging.Logger;
import java.util.zip.CRC32;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

/**
* Singleton helper class to add entries to zip archives.
* @since Feb 8, 2013
* @author edlunde
* 
*/
public class ZipAppendUtil {
	private static final Logger s_logger = Logger.getLogger( ZipAppendUtil.class.getPackage().getName() );
	private boolean busy;
	private static ZipAppendUtil s_instance;
	
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
	 * Extract contents of sourceArchive to a temporary file, append files to pathWithinArchive, and rezip.
	 * @param sourceArchive
	 * @param files
	 * @param pathWithinArchive
	 * @return
	 */
	public boolean addToZipArchive(File sourceArchive, File[] files, String pathWithinArchive){
		if (busy) return false;
		
		busy = true;
		if (files == null) return false;
		if(files.length >0){
			try{

				File tmpZip = File.createTempFile(sourceArchive.getName(), null);
				tmpZip.delete();
				if(!sourceArchive.renameTo(tmpZip)){
					s_logger.info("\nDictionary archive is occupied. Reading from a temporary copy instead.");
					copyFile(sourceArchive, tmpZip); //Since the chemcomp.zip file cannot be renamed while in use (Win only), copy it to tmpZip
				}
				byte[] buffer = new byte[1024];
				ZipInputStream zin = new ZipInputStream(new FileInputStream(tmpZip));
				ZipOutputStream out = new ZipOutputStream(new FileOutputStream(sourceArchive));
				for(int i = 0; i < files.length; i++){
					InputStream in = new FileInputStream(files[i]);
					ZipEntry entry =new ZipEntry(pathWithinArchive + files[i].getName());
					entry.setSize(files[i].length());
					entry.setTime(files[i].lastModified());
					CRC32 crc32 = new CRC32();
					out.putNextEntry(entry);
					for(int read = in.read(buffer); read > -1; read = in.read(buffer)){
						out.write(buffer, 0, read);
						crc32.update(buffer,0,read);
					}
					entry.setCrc(crc32.getValue());
					out.closeEntry();
					in.close();
				}
				for(ZipEntry ze = zin.getNextEntry(); ze != null; ze = zin.getNextEntry()){
					if(!zipEntryMatch(ze.getName(), files, pathWithinArchive)){
						ZipEntry destEntry = new ZipEntry(ze.getName());
						CRC32 crc32 = new CRC32();
						out.putNextEntry(destEntry);
						for(int read = zin.read(buffer); read > -1; read = zin.read(buffer)){
							out.write(buffer, 0, read);
							crc32.update(buffer,0,read);
						}

						destEntry.setCrc(crc32.getValue());
						out.closeEntry();
					}
				}
				out.close();
				zin.close();
				/*
				if(tmpZip.delete()) {
					s_logger.info("Successfully deleted temp zip file");
				}	
				s_logger.info("Wrote downloaded file(s) to "+sourceArchive.getAbsolutePath());
				*/
				busy = false;
				return true;
			}catch(Exception e){
				busy = false;
				s_logger.info(e.getMessage());
				return false;
			}
		}
		busy = false;
		return false;
	}
	
	/**
	 * Check for an entry.
	 * @param zeName
	 * @param files
	 * @param path
	 * @return
	 */
	private boolean zipEntryMatch(String zeName, File[] files, String path){
		for(int i = 0; i < files.length; i++){
			if((path + files[i].getName()).equals(zeName)){
				return true;
			}
		}
		return false;
	}

	/**
	 * Create a copy for a file.
	 * @param src
	 * @param dst
	 * @throws IOException
	 */
	public void copyFile(File src, File dst) throws IOException {
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
	        is.close();
	        os.close();
	    }
	}
	
	/**
	 * Check if the zip archive is being modified or referenced.
	 * @return
	 */
	public boolean isBusy() {
		return this.busy;
	}
	
	/* Prevent direct access to the constructor*/
	private ZipAppendUtil() {
		super();
	}
}
