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

import okhttp3.*;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.nio.channels.FileChannel;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

public class Download {

	//private static final Logger logger = LoggerFactory.getLogger(Download.class);
	private static final int DOWNLOAD_TIMEOUT_MS_DEFAULT = 60000;

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

        try (FileChannel source = new FileInputStream(src).getChannel(); FileChannel destination = new FileOutputStream(dst).getChannel()) {
            // we need the supress warnings here (the warning that the stream is not closed is harmless)
            // see http://stackoverflow.com/questions/12970407/does-filechannel-close-close-the-underlying-stream
            destination.transferFrom(source, 0, source.size());
        }
	}

	public static String getFileExtension(File f) {
		String fileName = f.getName();
		return fileName.substring(fileName.lastIndexOf('.') + 1);
	}

	public static String getFilePrefix(File f) {
		String fileName = f.getName();
		return fileName.substring(0, fileName.indexOf('.'));
	}

	static final OkHttpClient.Builder clients;
	static {
		/** https://square.github.io/okhttp/recipes/#response-caching-kt-java */
		int cacheSize = 512 * 1024 * 1024; //512m
		OkHttpClient.Builder bb = new OkHttpClient.Builder();
		bb.followRedirects(true);
		bb.followSslRedirects(true);
		bb.retryOnConnectionFailure(true);

//		try {
			File cacheDir = Paths.get(System.getProperty("java.io.tmpdir"), "biojava").toFile();
			cacheDir.mkdirs();
			Cache cache = new Cache(cacheDir, cacheSize);
			bb.cache(cache);

//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		clients = bb;
	}


//	/**
//	 * Converts path to Unix convention and adds a terminating slash if it was
//	 * omitted
//	 *
//	 * @param path original platform dependent path
//	 * @return path in Unix convention
//	 * @author Peter Rose
//	 * @since 3.2
//	 */
//	public static String toUnixPath(String path) {
//		String uPath = path;
//		if (uPath.contains("\\")) {
//			uPath = uPath.replaceAll("\\\\", "/");
//		}
//		// this should be removed, it's need since "\" is added AtomCache code
//		if (uPath.endsWith("//")) {
//			uPath = uPath.substring(0, uPath.length() - 1);
//		}
//		if (!uPath.endsWith("/")) {
//			uPath = uPath + "/";
//		}
//		return uPath;
//	}

	/**
	 * Expands ~ in paths to the user's home directory.
	 *
	 * <p>
	 * This does not work for some special cases for paths: Other users' homes
	 * (~user/...), and Tilde expansion within the path (/.../~/...)
	 *
	 * @param file
	 * @return
	 */
	public static String expandUserHome(String file) {
		if (file.startsWith("~" + File.separator)) {
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
	private static URLConnection prepareURLConnection(String url, int timeout) throws IOException {
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
	            if (e != null)
	                throw e;
	            Files.delete(dir);
	            return FileVisitResult.CONTINUE;
	        }
	    });
	}

	public static InputStream stream(URL url) throws IOException {
		return stream(url, DOWNLOAD_TIMEOUT_MS_DEFAULT);
	}

	/**
	 * Open a URL and return an InputStream to it
	 * if acceptGzipEncoding == true, use GZIPEncoding to
	 * compress communication.
	 * <p>
	 * The caller is responsible to close the returned InputStream not to cause
	 * resource leaks.
	 * @param url the input URL to be read
	 * @param timeoutMS
	 * @return an {@link InputStream} of response
	 * @throws IOException due to an error opening the URL
	 */
	public static InputStream stream(URL url, int timeoutMS) throws IOException {
//		InputStream inStream;
//		URLConnection huc = URLConnectionTools.openURLConnection(url,timeout);
//
//		if ( acceptGzipEncoding) huc.setRequestProperty("Accept-Encoding", "gzip");
//
//		String contentEncoding = huc.getContentEncoding();
//
//		inStream = huc.getInputStream();

		OkHttpClient client = clients.build(); //TODO ThreadLocal

		Call c = client.newCall(new Request.Builder().url(url)
				//.cacheControl(CacheControl.FORCE_CACHE)
				.build());

		c.timeout().deadline(timeoutMS, TimeUnit.MILLISECONDS);

		try (Response r = c.execute()) {
			InputStream inStream;
			ResponseBody b = r.body();
			if (!r.isSuccessful()) {
				//throw new IOException("Unexpected code " + r);
				inStream = new ByteArrayInputStream(new byte[] { } ); //HACK
			} else {
				inStream =
						//b.byteStream();
						new ByteArrayInputStream(b.bytes()); //b.byteStream();
			}

			if (b.contentType().toString().contains("-gzip"))
				return new GZIPInputStream(inStream);

			return inStream;
		}


	}

	/**
	 * Download the content provided at URL url and store the result to a local
	 * file, using a temp file to cache the content in case something goes wrong
	 * in download
	 *
	 * @param url
	 * @param destination
	 * @throws IOException
	 */
	public static void downloadFile(URL url, File destination) throws IOException {


		OkHttpClient client = clients.build(); //TODO ThreadLocal

		try (Response r = client.newCall(new Request.Builder().url(url).build()).execute()) {
			if (!r.isSuccessful()) {
				//throw new IOException("Unexpected code " + r);
				return;
			}

			ResponseBody b = r.body();

			FileChannel dst = new FileOutputStream(destination).getChannel();
			dst.transferFrom(b.source(), 0, b.contentLength());

			System.out.println("Response 1 response:          " + r);
			System.out.println("Response 1 cache response:    " + r.cacheResponse());
			System.out.println("Response 1 network response:  " + r.networkResponse());
		}


//		int count = 0;
//		int maxTries = 10;
//		int timeout = 60000; //60 sec
//
//		File tempFile = File.createTempFile(getFilePrefix(destination), "." + getFileExtension(destination));
//
//		// Took following recipe from stackoverflow:
//		// http://stackoverflow.com/questions/921262/how-to-download-and-save-a-file-from-internet-using-java
//		// It seems to be the most efficient way to transfer a file
//		// See: http://docs.oracle.com/javase/7/docs/api/java/nio/channels/FileChannel.html
//		ReadableByteChannel rbc = null;
//		FileOutputStream fos = null;
//		while (true) {
//			try {
//				URLConnection connection = prepareURLConnection(url.toString(), timeout);
//				connection.connect();
//
//				rbc = Channels.newChannel(connection.getInputStream());
//				fos = new FileOutputStream(tempFile);
//				fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
//				break;
//			} catch (SocketTimeoutException e) {
//				if (++count == maxTries) throw e;
//			} finally {
//				if (rbc != null) {
//					rbc.close();
//				}
//				if (fos != null) {
//					fos.close();
//				}
//			}
//		}
//
//		logger.debug("Copying temp file {} to final location {}", tempFile, destination);
//		copy(tempFile, destination);
//
//		// delete the tmp file
//		tempFile.delete();

	}

//	/**
//	 * Recursively delete a folder & contents
//	 *
//	 * @param dir directory to delete
//	 */
//	public static void deleteDirectory(String dir) throws IOException {
//		deleteDirectory(Paths.get(dir));
//	}
//
//
//	public static void main(String[] args) {
//		String url;
//		url = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/";
//		System.out.format("%s\t%s%n", ping(url, 200), url);
//		url = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/foo";
//		System.out.format("%s\t%s%n", ping(url, 200), url);
//		url = "http://scopzzz.mrc-lmb.cam.ac.uk/scop/parse/";
//		System.out.format("%s\t%s%n", ping(url, 200), url);
//		url = "scop.mrc-lmb.cam.ac.uk";
//		System.out.format("%s\t%s%n", ping(url, 200), url);
//		url = "http://scop.mrc-lmb.cam.ac.uk";
//		System.out.format("%s\t%s%n", ping(url, 200), url);
//	}

}
