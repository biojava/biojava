/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These
 * should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the
 * biojava-l mailing list, visit the home page at:
 *
 * http://www.biojava.org/
 *
 * Created on Feb 23, 2012 Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.core.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.HttpURLConnection;
import java.net.SocketTimeoutException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Scanner;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FileDownloadUtils {

	private static final String SIZE_EXT = ".size";
	private static final String HASH_EXT = ".hash";
	private static final Logger logger = LoggerFactory.getLogger(FileDownloadUtils.class);

	public enum Hash{
		MD5, SHA1, SHA256, UNKNOWN
	}

	/**
	 * Gets the file extension of a file, excluding '.'.
	 * If the file name has no extension the file name is returned.
	 * @param f a File
	 * @return The extension
	 */
	public static String getFileExtension(File f) {
		String fileName = f.getName();
		String ext = "";
		int mid = fileName.lastIndexOf(".");
		ext = fileName.substring(mid + 1);
		return ext;
	}

	/**
	 * Gets the file name up to and excluding the first
	 * '.' character. If there is no extension, the full filename
	 * is returned.
	 * @param f A file
	 * @return A possibly empty but non-null String.
	 */
	public static String getFilePrefix(File f) {
		String fileName = f.getName();
		int mid = fileName.indexOf(".");
		if (mid < 0) {
			return fileName;
		}
		return fileName.substring(0, mid);
	}

	/**
	 * Download the content provided at URL url and store the result to a local
	 * file, using a temp file to cache the content in case something goes wrong
	 * in download. A timeout of 60 seconds is hard-coded and 10 retries are attempted.
	 *
	 * @param url
	 * @param destination
	 * @throws IOException
	 */
	public static void downloadFile(URL url, File destination) throws IOException {
		int count = 0;
		int maxTries = 10;
		int timeout = 60000; //60 sec

		File tempFile = Files.createTempFile(getFilePrefix(destination), "." + getFileExtension(destination)).toFile();

		// Took following recipe from stackoverflow:
		// http://stackoverflow.com/questions/921262/how-to-download-and-save-a-file-from-internet-using-java
		// It seems to be the most efficient way to transfer a file
		// See: http://docs.oracle.com/javase/7/docs/api/java/nio/channels/FileChannel.html
		ReadableByteChannel rbc = null;
		FileOutputStream fos = null;
		while (true) {
			try {
				URLConnection connection = prepareURLConnection(url.toString(), timeout);
				connection.connect();
				InputStream inputStream = connection.getInputStream();

				rbc = Channels.newChannel(inputStream);
				fos = new FileOutputStream(tempFile);
				fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
				break;
			} catch (SocketTimeoutException e) {
				if (++count == maxTries) throw e;
			} finally {
				if (rbc != null) {
					rbc.close();
				}
				if (fos != null) {
					fos.close();
				}
			}
		}

		logger.debug("Copying temp file [{}] to final location [{}]", tempFile, destination);
		Files.copy(tempFile.toPath(), destination.toPath(), StandardCopyOption.REPLACE_EXISTING);

		// delete the tmp file
		tempFile.delete();

	}
	
	/**
	 * Creates validation files beside a file to be downloaded.<br>
	 * Whenever possible, for a <code>file.ext</code> file, it creates 
	 * <code>file.ext.size</code> and <code>file.hash</code> for in the same 
	 * folder where <code>file.ext</code> exists.
	 * If the file connection size could not be deduced from the URL, no size file is created. 
	 * If <code>hashURL</code> is <code>null</code>, no hash file is created.
	 * @param url the remote file URL to download
	 * @param localDestination the local file to download into
	 * @param hashURL the URL of the hash file to download. Can be <code>null</code>.
	 * @param hash The Hashing algorithm. Ignored if <code>hashURL</code> is <code>null</code>.
	 */
	public static void createValidationFiles(URL url, File localDestination, URL hashURL, Hash hash){
		try {
			URLConnection resourceConnection = url.openConnection();
			createValidationFiles(resourceConnection, localDestination, hashURL, FileDownloadUtils.Hash.UNKNOWN);
		} catch (IOException e) {
			logger.warn("could not open connection to resource file due to exception: {}", e.getMessage());
		}
	}
	/**
	 * Creates validation files beside a file to be downloaded.<br>
	 * Whenever possible, for a <code>file.ext</code> file, it creates 
	 * <code>file.ext.size</code> and <code>file.hash_XXXX</code> in the same 
	 * folder where <code>file.ext</code> exists (XXXX may be DM5, SHA1, or SHA256).
	 * If the file connection size could not be deduced from the resourceUrlConnection 
	 * {@link URLConnection}, no size file is created. 
	 * If <code>hashURL</code> is <code>null</code>, no hash file is created.<br>
	 * <b>N.B.</b> None of the hashing algorithms is implemented (yet), because we did not need any of them yet.
	 * @param resourceUrlConnection the remote file URLConnection to download
	 * @param localDestination the local file to download into
	 * @param hashURL the URL of the hash file to download. Can be <code>null</code>.
	 * @param hash The Hashing algorithm. Ignored if <code>hashURL</code> is <code>null</code>.
	 * @since 7.0.0
	 */
	public static void createValidationFiles(URLConnection resourceUrlConnection, File localDestination, URL hashURL, Hash hash){
		long size = resourceUrlConnection.getContentLengthLong();
		if(size == -1) {
			logger.debug("Could not find expected file size for resource {}. Size validation metadata file won't be available for this download.", resourceUrlConnection.getURL());
		} else {
			logger.debug("Content-Length: {}", size);
			File sizeFile = new File(localDestination.getParentFile(), localDestination.getName() + SIZE_EXT);
			try (PrintStream sizePrintStream = new PrintStream(sizeFile)) {
				sizePrintStream.print(size);
			} catch (FileNotFoundException e) {
				logger.warn("Could not write size validation metadata file due to exception: {}", e.getMessage());
			}
		}
		
		if(hashURL == null)
			return;

		if(hash == Hash.UNKNOWN)
			throw new IllegalArgumentException("Hash URL given but algorithm is unknown");
		try {
			File hashFile = new File(localDestination.getParentFile(), String.format("%s%s_%s", localDestination.getName(), HASH_EXT, hash));
			downloadFile(hashURL, hashFile);
		} catch (IOException e) {
			logger.warn("Could not write validation hash file due to exception: {}", e.getMessage());
		}
	}
	
	/**
	 * Validate a local file based on pre-existing metadata files for size and hash.<br>
	 * If the passed in <code>localFile</code> parameter is a file named <code>file.ext</code>, the function searches in the same folder for:
	 * <ul>
	 * <li><code>file.ext.size</code>: If found, it compares the size stored in it to the length of <code>localFile</code> (in bytes).</li>
	 * <li><code>file.ext.hash_XXXX (where XXXX is DM5, SHA1, or SHA256)</code>: If found, it compares the size stored in it to the hash code of <code>localFile</code>.</li>
	 * </ul>
	 * If any of these comparisons fail, the function returns <code>false</code>. otherwise it returns true.
	 * <p>
	 * <b>N.B.</b> None of the 3 common verification hashing algorithms are implement yet.
	 * @param localFile The file to validate
	 * @return <code>false</code> if any of the size or hash code metadata files exists but its contents does not match the expected value in the file, <code>true</code> otherwise.
	 * @since 7.0.0
	 */
	public static boolean validateFile(File localFile) {
		File sizeFile = new File(localFile.getParentFile(), localFile.getName() + SIZE_EXT);
		if(sizeFile.exists()) {
            try (Scanner scanner = new Scanner(sizeFile)) {
                long expectedSize = scanner.nextLong();
                long actualSize = localFile.length();
                if (expectedSize != actualSize) {
                    logger.warn("File [{}] size ({}) does not match expected size ({}).", localFile, actualSize, expectedSize);
                    return false;
                }
            } catch (FileNotFoundException e) {
                logger.warn("could not validate size of file [{}] because no size metadata file exists.", localFile);
            }
		}

		File[] hashFiles = localFile.getParentFile().listFiles(new FilenameFilter() {
			final String hashPattern = String.format("%s%s_(%s|%s|%s)", localFile.getName(), HASH_EXT, Hash.MD5, Hash.SHA1, Hash.SHA256);
			@Override
			public boolean accept(File dir, String name) {
				return name.matches(hashPattern);
			}
		});
		if(hashFiles.length > 0) {
			File hashFile = hashFiles[0];
			String name = hashFile.getName();
			String algo = name.substring(name.lastIndexOf('_') + 1);
			switch (Hash.valueOf(algo)) {
			case MD5:
			case SHA1:
			case SHA256:
				throw new UnsupportedOperationException("Not yet implemented");
			case UNKNOWN:
			default: // No need. Already checked above
				throw new IllegalArgumentException("Hashing algorithm not known: " + algo);
			}
		}
		
		return true;
	}

	/**
	 * Converts path to Unix convention and adds a terminating slash if it was
	 * omitted. 
	 *
	 * @param path original platform dependent path
	 * @return path in Unix convention
	 * @author Peter Rose
	 * @since 3.2
	 */
	public static String toUnixPath(String path) {
		String uPath = path;
		if (uPath.contains("\\")) {
			uPath = uPath.replaceAll("\\\\", "/");
		}
		// this should be removed, it's need since "\" is added AtomCache code
		if (uPath.endsWith("//")) {
			uPath = uPath.substring(0, uPath.length() - 1);
		}
		if (!uPath.endsWith("/")) {
			uPath = uPath + "/";
		}
		return uPath;
	}

	/**
	 * Expands ~ in paths to the user's home directory.
	 *
	 * <p>
	 * This does not work for some special cases for paths: Other users' homes
	 * (~user/...), and Tilde expansion within the path (/.../~/...). In these cases
	 *  the original argument is returned.
	 *
	 * @param file A filepath starting with a tilde
	 * @return An absolute path
	 */
	public static String expandUserHome(String file) {
		// replace any / with the proper separator (/ or \ for Linux and Windows respectively).
		file = file.replaceAll("/", "\\"+File.separator); //The "\\" is to escape the separator if needed.
		if (file.startsWith("~") && (file.length() == 1 || File.separator.equals(file.substring(1, 2)))) {
			file = System.getProperty("user.home") + file.substring(1);
		}
		return file;
	}

	/**
	 * Pings a HTTP URL. This effectively sends a HEAD request and returns
	 * <code>true</code> if the response code is in the 200-399 range.
	 *
	 * @param url The HTTP URL to be pinged.
	 * @param timeout The timeout in millis for both the connection timeout and
	 * the response read timeout. Note that the total timeout is effectively two
	 * times the given timeout.
	 * @return <code>true</code> if the given HTTP URL has returned response
	 * code 200-399 on a HEAD request within the given timeout, otherwise
	 * <code>false</code>.
	 * @author BalusC,
	 * http://stackoverflow.com/questions/3584210/preferred-java-way-to-ping-a-http-url-for-availability
	 */
	public static boolean ping(String url, int timeout) {
		//url = url.replaceFirst("https", "http"); // Otherwise an exception may be thrown on invalid SSL certificates.

		try {
			HttpURLConnection connection = (HttpURLConnection) prepareURLConnection(url, timeout);
			connection.setRequestMethod("HEAD");
			int responseCode = connection.getResponseCode();
			return (200 <= responseCode && responseCode <= 399);
		} catch (IOException exception) {
			return false;
		}
	}

	/**
	 * Prepare {@link URLConnection} with customised timeouts.
	 *
	 * @param url The URL
	 * @param timeout The timeout in millis for both the connection timeout and
	 * the response read timeout. Note that the total timeout is effectively two
	 * times the given timeout.
	 *
	 * <p>
	 * Example of code.      <code>
		 * UrlConnection conn = prepareURLConnection("http://www.google.com/", 20000);
	 * conn.connect();
	 * conn.getInputStream();
	 * </code>
	 * <p>
	 *
	 * <strong>NB. User should execute connect() method before getting input
	 * stream.</strong>
	 * @return
	 * @throws IOException
	 * @author Jacek Grzebyta
	 */
	public static URLConnection prepareURLConnection(String url, int timeout) throws IOException {
		URLConnection connection = new URL(url).openConnection();
		connection.setReadTimeout(timeout);
		connection.setConnectTimeout(timeout);
		return connection;
	}

	/**
	 * Recursively delete a folder &amp; contents
	 *
	 * @param dir directory to delete
	 */
	public static void deleteDirectory(Path dir) throws IOException {
		if(dir == null || !Files.exists(dir))
			return;
		Files.walkFileTree(dir, new SimpleFileVisitor<>() {
	        @Override
	        public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
	            Files.delete(file);
	            return FileVisitResult.CONTINUE;
	        }

	        @Override
	        public FileVisitResult postVisitDirectory(Path dir, IOException e) throws IOException {
	            if (e != null) {
	                throw e;
	            }
	            Files.delete(dir);
	            return FileVisitResult.CONTINUE;
	        }
	    });
	}
	/**
	 * Recursively delete a folder &amp; contents
	 *
	 * @param dir directory to delete
	 */
	public static void deleteDirectory(String dir) throws IOException {
		deleteDirectory(Paths.get(dir));
	}

}
