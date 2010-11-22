package org.biojava.bio.structure.align.util;



import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class SynchronizedOutFile {

	File file;

	String[] tmp;

	int ARR_SIZE=100;
	int counter;

	/** create a thread safe wrapper for working with this file
	 * 
	 * @param f
	 */
	public SynchronizedOutFile(File f) throws FileNotFoundException, IOException{
		if ( f.isDirectory())
			throw new FileNotFoundException("please provide a file and not a directory");

		if ( ! f.exists()){
			f.createNewFile();
		}
		file = f;
		tmp = new String[ARR_SIZE];
		counter = -1;
	}

	public synchronized void write(String message) throws IOException{ 
		counter++;
		tmp[counter] = message;
		if (counter >= ARR_SIZE - 1 ) {
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
		try {
			out = new BufferedOutputStream(new  FileOutputStream(file, true));
			for ( int i = 0 ; i <= counter ; i++){
				byte data[] = tmp[i].getBytes();
				out.write(data, 0, data.length);
			}

		} catch (IOException x) {
			System.err.println(x);
		} finally {
			if (out != null) {
				out.flush();
				out.close();
			}
		}
	}


}
