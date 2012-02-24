/**
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
 * Created on Feb 23, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.io.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;

import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class FileDownloadUtils {

	/** Copy the content of file A to B
	 * 
	 * @param src
	 * @param dst
	 * @throws IOException
	 */
	public static void copy(File src, File dst) throws IOException {

		InputStream in = new FileInputStream(src);
		OutputStream out = new FileOutputStream(dst);

		// Transfer bytes from in to out
		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0) {
			out.write(buf, 0, len);
		}
		in.close();
		out.close();
	}

	public static String getFileExtension(File f){
		String fileName = f.getName();
		//String fname="";
		String ext="";
		int mid= fileName.lastIndexOf(".");
		//fname=fileName.substring(0,mid);
		ext=fileName.substring(mid+1,fileName.length());  
		//System.out.println("File name ="+fname);
		//System.out.println("Extension ="+ext);
		return ext;
	}

	public static String getFilePrefix(File f){
		String fileName = f.getName();
		String fname="";

		int mid= fileName.indexOf(".");
		fname=fileName.substring(0,mid);

		return fname;
	}


	/** Download the content provided at URL url and stores the result to a local file
	 * 
	 * @param url
	 * @param destination
	 * @throws IOException
	 */
	public static void downloadGzipCompressedFile(URL url, File destination) throws IOException{


		InputStream uStream = url.openStream();
		InputStream conn = new GZIPInputStream(uStream);

		File tempFile  = File.createTempFile(getFilePrefix(destination), "."+ getFileExtension(destination));

		try {
			System.out.println("downloading " + url + " to " + tempFile.getAbsolutePath());
			FileOutputStream outPut = new FileOutputStream(tempFile);
			GZIPOutputStream gzOutPut = new GZIPOutputStream(outPut);
			PrintWriter pw = new PrintWriter(gzOutPut);

			BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(conn));
			String line;
			while ((line = fileBuffer.readLine()) != null) {
				pw.println(line);
			}
			pw.flush();
			pw.close();

			outPut.flush();
			outPut.close();
			conn.close();
			uStream.close();
		} catch (Exception e) {
			e.printStackTrace();
			if ( conn != null)
				conn.close();
			if (uStream != null) 
				uStream.close();
			throw new IOException(e.getMessage());
			
		}
		// copy file name to **real** location (without the tmpFileName)
		// prepare destination
		System.out.println("copying to " + destination);

		copy(tempFile, destination);

		// delete the tmp file			
		tempFile.delete();

	}

	public static File downloadFileIfAvailable(URL url, File destination) throws IOException {

		InputStream uStream = null;
		InputStream conn = null;
		try {
			uStream = url.openStream();
			conn = new GZIPInputStream(uStream);
		} catch (IOException e1) {
			System.err.println("Problem while downloading file " + url  );
			//e1.printStackTrace();
			try {
				if (uStream != null) {
					uStream.close();
				}	

				if (conn != null) {
					conn.close();
				}

			} catch (IOException e) {
				e.printStackTrace();
			}
			return null;
		} 


		FileOutputStream outPut = null;
		GZIPOutputStream gzOutPut = null;
		File tempFile  = File.createTempFile(getFilePrefix(destination), "."+ getFileExtension(destination));
		try {


			outPut = new FileOutputStream(tempFile);
			gzOutPut = new GZIPOutputStream(outPut);
			PrintWriter pw = new PrintWriter(gzOutPut);

			BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(conn));
			String line;
			while ((line = fileBuffer.readLine()) != null) {
				pw.println(line);
			}
			pw.flush();
			pw.close();

			outPut.flush();
			outPut.close();
			conn.close();
			uStream.close();

		} catch (Exception e){
			System.err.println("Problem while downloading " + url );
			e.printStackTrace();
			return null;
		} finally {	
			if ( conn != null){
				try {

					conn.close();
				} catch (IOException e){
					e.printStackTrace();
				}
			}
			if ( uStream != null){
				try { 
					uStream.close();					
				}catch (IOException e){
					e.printStackTrace();
				}
			}
			try {
				if (outPut != null) {
					outPut.close();
				}
				if (gzOutPut != null) {
					gzOutPut.close();
				}
			} catch (IOException e) {
			}
		}
		System.out.println("Writing to " + destination);

		copy(tempFile, destination);

		return destination;
	}

	/**
	 * Converts path to Unix convention and adds a terminating slash if it was omitted
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
			uPath = uPath.substring(0, uPath.length()-1);
		}
		if (! uPath.endsWith("/")) {
			uPath = uPath + "/";
		}
		return uPath;
	}



}
