package org.biojava.bio.structure.server;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.MMCIFFileReader;

import org.biojava.bio.structure.io.mmcif.MMcifParser;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava3.core.util.InputStreamProvider;

/**
 * 
 * @deprecated
 *
 */
public class MMCIFFileInstallation implements PDBInstallation
{

	public static final Logger logger = Logger.getLogger("org.biojava.bio.structure");
	private File filePath;
	List<String> extensions  ;
	MMCIFFileReader reader ;
	boolean autoFetch;

	public MMCIFFileInstallation(File filePath){

		if (! filePath.isDirectory()){
			throw new IllegalArgumentException("the provided path does not point to a directory!");
		}

		this.filePath = filePath;
		reader = new MMCIFFileReader();
		extensions    = new ArrayList<String>();

		extensions.add(".cif");
		extensions.add(".mmcif");
		extensions.add(".cif.gz");
		extensions.add(".mmcif.gz");
		extensions.add(".cif.Z");
		extensions.add(".mmcif.Z");

		autoFetch = false;
	}

	/** should the parser to fetch missing mmCif files from the RCSB FTP server automatically?
	 *  default is false
	 * @return flag
	 */
	public boolean isAutoFetch() {
		return autoFetch;
	}

	/** tell the parser to fetch missing mmCif files from the RCSB FTP server automatically.
	 *
	 * default is false. If true, new PDB files will be automatically stored in the Path and gzip compressed.
	 *
	 * @param autoFetch
	 */
	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}


	private InputStream getInputStream(String pdbId)
	throws IOException
	{
		InputStream inputStream =null;

		String pdbFile = null ;
		File f = null ;

		// this are the possible PDB file names...
		String fpath = filePath+"/"+pdbId;
		String ppath = filePath +"/pdb"+pdbId;

		String split = pdbId.substring(1,3);
		String fsplit = filePath+"/"+split+"/"+pdbId;

		String[] paths = new String[]{fpath,ppath,fsplit};
		for (String testpath : paths) {
			for (int i=0 ; i<extensions.size();i++){
				String ex = (String)extensions.get(i) ;
				//System.out.println("PDBFileReader testing: "+testpath+ex);
				f = new File(testpath+ex) ;

				if ( f.exists()) {
					//System.out.println("found!");
					pdbFile = testpath+ex ;

					InputStreamProvider isp = new InputStreamProvider();

					inputStream = isp.getInputStream(pdbFile);
					break;
				}

				if ( pdbFile != null) break;
			}
		}

		if ( pdbFile == null ) {
			if ( autoFetch)
				return downloadAndGetInputStream(pdbId);

			String message = "no structure with PDB code " + pdbId + " found!" ;
			throw new IOException (message);
		}

		return inputStream ;

	}

	private InputStream downloadAndGetInputStream(String pdbId)
	throws IOException{
		//PDBURLReader reader = new PDBURLReader();
		//Structure s = reader.getStructureById(pdbId);
		File tmp = downloadCIF(pdbId);
		if ( tmp != null ) {
			InputStreamProvider prov = new InputStreamProvider();
			return prov.getInputStream(tmp);


		} else {
			throw new IOException("could not find PDB " + pdbId + " in file system and also could not download");
		}
	}

	public File downloadCIF(String pdbId){
		File tempFile = new File(filePath+"/"+pdbId+".cif.gz");
		File pdbHome  = new File(filePath.getAbsolutePath());

		if ( ! pdbHome.canWrite() ){
			System.err.println("can not write to " + pdbHome);
			return null;
		}

		String ftp = String.format("ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/mmCIF/%s.cif.gz", pdbId.toLowerCase());

		System.out.println("Fetching " + ftp);
		try {
			URL url = new URL(ftp);
			InputStream conn = url.openStream();

			// prepare destination
			System.out.println("writing to " + tempFile);


			FileOutputStream outPut = new FileOutputStream(tempFile);
			GZIPOutputStream gzOutPut = new GZIPOutputStream(outPut);
			PrintWriter pw = new PrintWriter(gzOutPut);

			BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(new GZIPInputStream(conn)));
			String line;
			while ((line = fileBuffer.readLine()) != null) {
				pw.println(line);
			}
			pw.flush();
			pw.close();
			outPut.close();
			conn.close();
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return tempFile;
	}

	public void addPDBFilter(PDBFilter filter) {
		// TODO Auto-generated method stub

	}

	public void clearFilters() {
		// TODO Auto-generated method stub

	}

	public List<PDBHeader> getAll() {
		// TODO Auto-generated method stub
		return null;
	}

	public PDBHeader getPDBHeader(String pdbId) {
		// TODO Auto-generated method stub
		return null;
	}

	public Structure getStructure(String pdbId)  {

		Structure struc = null;


		try {
			InputStream inStream = getInputStream(pdbId);
			BufferedReader buf ;
			try {
				buf = getBufferedReader(inStream);

			} catch (IOException e) {
				e.printStackTrace();
				throw new IOException ("error initializing BufferedReader");
			}
			//System.out.println("done");


			MMcifParser pdbpars = new SimpleMMcifParser();
			SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
			pdbpars.addMMcifConsumer(consumer);
			pdbpars.parse(buf) ;
			struc = consumer.getStructure();
		} catch(IOException e){
			e.printStackTrace();
		}
		return struc ;


	}

	private BufferedReader getBufferedReader(InputStream inStream)
	throws IOException {

		BufferedReader buf ;
		if (inStream == null) {
			throw new IOException ("input stream is null!");
		}

		buf = new BufferedReader (new InputStreamReader (inStream));
		return buf ;

	}

	public boolean hasNext() {
		// TODO Auto-generated method stub
		return false;
	}

	public Structure next() {
		// TODO Auto-generated method stub
		return null;
	}

}
