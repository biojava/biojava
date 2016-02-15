package org.biojava.nbio.structure.rcsb;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RCSBUpdates {
	
	// The URL for acquiring the data
	public static final String baseURL = "ftp://ftp.rcsb.org/pub/pdb/data/status/latest/";

	/**
	 * 
	 * @return A map mapping each field (defined by a separate FTP file) to the PDB ids in the field. The possible fields
	 * are: added.models, added.nmr, added.pdb, added.sf, modified.cs, modified.models, modified.nmr, modified.pdb, modified.sf, 
	 * obsolete.cs, obsolete.models, obsolete.nmr, obsolete.pdb, obsolete.sf
	 * @throws IOException 
	 */
	public Map<String,String[]> getUpdates() throws IOException{
		
		Map<String,String[]> outMap = new HashMap<String, String[]>();
		// A list of files to get
		String[] newStringList = {"added.models","added.nmr","added.pdb","added.sf","modified.cs","modified.models",
				"modified.nmr","modified.pdb","modified.sf","obsolete.cs","obsolete.models","obsolete.nmr","obsolete.pdb","obsolete.sf"};
		for(String fileName: newStringList){
			String[] thisList = readURL(baseURL+""+fileName);
			outMap.put(fileName, thisList);
		}		
		return outMap;

	}


	/**
	 * 
	 * @param urlIn The url to be read
	 * @return A list of PDB ids as strings
	 * @throws IOException
	 */
	private String[] readURL(String urlIn) throws IOException{
		List<String> outList = new ArrayList<String>();
		// create a url object
		URL url = new URL(urlIn);

		// create a urlconnection object
		URLConnection urlConnection = url.openConnection();

		// wrap the urlconnection in a bufferedreader
		try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(urlConnection.getInputStream()))) {

			String line;

			// read from the urlconnection via the bufferedreader
			while ((line = bufferedReader.readLine()) != null)
			{
				outList.add(line);
			}

		}

		return outList.toArray(new String[outList.size()]);
	}
}
