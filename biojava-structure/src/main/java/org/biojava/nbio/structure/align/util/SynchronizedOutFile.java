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
package org.biojava.nbio.structure.align.util;


import java.io.*;
import java.util.zip.GZIPOutputStream;

public class SynchronizedOutFile {

	File file;

	String[] tmp;

	int ARR_SIZE=100;
	Integer counter;

	boolean useGzipCompression = false;


	/** Create a thread safe wrapper for writing to this file, the file will be gzip compressed.
	 *
	 * @param f file to write to
	 * @param gzipCompress flag if file should be gzip compressed
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public SynchronizedOutFile(File f, boolean gzipCompress) throws FileNotFoundException, IOException{
		if ( f.isDirectory())
			throw new FileNotFoundException("please provide a file and not a directory");

		if ( ! f.exists()){
			System.out.println("creating output file: " + f.getAbsolutePath());
			f.createNewFile();
		}
		file = f;
		tmp = new String[ARR_SIZE];
		counter = -1;
		useGzipCompression = gzipCompress;

	}

	/** create a thread safe wrapper for working with this file
	 *
	 * @param f
	 */
	public SynchronizedOutFile(File f) throws FileNotFoundException, IOException{

		this(f,false);

	}

	public synchronized void write(String message) throws IOException{

		synchronized (counter){
			counter++;
			tmp[counter] = message;
			if (counter >= ARR_SIZE - 1 ) {
				writeArr();
				counter = -1;
			}
		}


	}

	public synchronized void flush() throws IOException {
		synchronized (counter){
			writeArr();
			counter = -1;
		}
	}

	public void close() throws IOException{
		writeArr();
		tmp = new String[ARR_SIZE];
	}

	private void writeArr() throws IOException{


		OutputStream out = null;
		FileOutputStream fileOutputStream=null;
		try {
			//This is less code-redundant
			fileOutputStream = new FileOutputStream(file, true);
			OutputStream outputstream = useGzipCompression? new GZIPOutputStream(fileOutputStream) : fileOutputStream;
			out = new BufferedOutputStream(outputstream);

			for ( int i = 0 ; i <= counter ; i++){
				if ( tmp[i] == null )
					continue;
				byte[] data = tmp[i].getBytes();
				out.write(data, 0, data.length);
			}

		} catch (Exception x) {
			System.err.println(x);
		} finally {
			if (out != null) {
				out.flush();
				out.close();
			}
		}
	}


}
