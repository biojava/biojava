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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.SocketTimeoutException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.channels.ReadableByteChannel;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FileDownloadUtils {

	private static final Logger logger = LoggerFactory.getLogger(FileDownloadUtils.class);

	/**
	 * Copy the content of file src to dst TODO since java 1.7 this is provided
	 * in java.nio.file.Files
	 *
	 * @param src
	 * @param dst
	 * @throws IOException
	 */
	@SuppressWarnings("resource")
	public static void copy(File src, File dst) throws IOException {

		// Took following recipe from
		// http://stackoverflow.com/questions/106770/standard-concise-way-to-copy-a-file-in-java
		// The nio package seems to be the most efficient way to copy a file
		FileChannel source = null;
		FileChannel destination = null;

		try {
			// we need the supress warnings here (the warning that the stream is not closed is harmless)
			// see http://stackoverflow.com/questions/12970407/does-filechannel-close-close-the-underlying-stream
			source = new FileInputStream(src).getChannel();
			destination = new FileOutputStream(dst).getChannel();
			destination.transferFrom(source, 0, source.size());
		} finally {
			if (source != null) {
				source.close();
			}
			if (destination != null) {
				destination.close();
			}
		}
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
		ext = fileName.substring(mid + 1, fileName.length());
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

		File tempFile = File.createTempFile(getFilePrefix(destination), "." + getFileExtension(destination));

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

		logger.debug("Copying temp file {} to final location {}", tempFile, destination);
		copy(tempFile, destination);

		// delete the tmp file
		tempFile.delete();

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
	 * <bold>NB. User should execute connect() method before getting input
	 * stream.</bold>
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
	 * Recursively delete a folder & contents
	 *
	 * @param dir directory to delete
	 */
	public static void deleteDirectory(Path dir) throws IOException {
		if(dir == null || !Files.exists(dir))
			return;
		Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
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
	 * Recursively delete a folder & contents
	 *
	 * @param dir directory to delete
	 */
	public static void deleteDirectory(String dir) throws IOException {
		deleteDirectory(Paths.get(dir));
	}

}