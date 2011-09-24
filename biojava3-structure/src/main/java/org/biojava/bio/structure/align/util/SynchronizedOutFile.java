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
	Integer counter;

	/** create a thread safe wrapper for working with this file
	 * 
	 * @param f
	 */
	public SynchronizedOutFile(File f) throws FileNotFoundException, IOException{
		if ( f.isDirectory())
			throw new FileNotFoundException("please provide a file and not a directory");

		if ( ! f.exists()){
			System.out.println("creating output file: " + f.getAbsolutePath());
			f.createNewFile();
		}
		file = f;
		tmp = new String[ARR_SIZE];
		counter = -1;
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

	public void close() throws IOException{
		writeArr();
		tmp = new String[ARR_SIZE];
	}

	private void writeArr() throws IOException{

		
		OutputStream out = null;
		try {
			out = new BufferedOutputStream(new  FileOutputStream(file, true));
			for ( int i = 0 ; i <= counter ; i++){
				if ( tmp[i] == null )
					continue;
				byte data[] = tmp[i].getBytes();
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
