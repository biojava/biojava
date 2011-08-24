package org.biojava.bio.structure.io;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava3.core.util.InputStreamProvider;


/** The "Sandbox" style of organizing files is  to have a directory structure like below, i.e. the files are organized into
 * <ul>
 * 	<li>directory with two characters, based on the two middle characters of a PDB ID</li>
 * 	<li>directory of PDB ID</li>
 * 	<li>several files that are available for this PDB ID</li>
 * </ul>
 * 
 * <pre>
a1/2a1v/2a1v.cif.gz
a1/2a1v/2a1v.dssp.gz
a1/2a1v/2a1v.pdb-250.jpg.gz
a1/2a1v/2a1v.pdb-500.jpg.gz
a1/2a1v/2a1v.pdb-65.jpg.gz
a1/2a1v/2a1v.pdb-80.jpg.gz
a1/2a1v/2a1v.pdb1-250.jpg.gz
a1/2a1v/2a1v.pdb1-500.jpg.gz
a1/2a1v/2a1v.pdb1-65.jpg.gz
a1/2a1v/2a1v.pdb1-80.jpg.gz
a1/2a1v/2a1v.pdb1.gz
a1/2a1v/2a1v.stride.gz
a1/2a1v/2a1v.xml.gz
a1/2a1v/pdb2a1v.ent.gz
a1/2a1v/r2a1vsf.ent.gz
a1/2a1w/2a1w-deriv.cif.gz
a1/2a1w/2a1w-extatom.xml.gz
a1/2a1w/2a1w-noatom.xml.gz
a1/2a1w/2a1w.cif.gz
a1/2a1w/2a1w.dssp.gz
a1/2a1w/2a1w.pdb-250.jpg.gz
a1/2a1w/2a1w.pdb-500.jpg.gz
a1/2a1w/2a1w.pdb-65.jpg.gz
a1/2a1w/2a1w.pdb-80.jpg.gz
a1/2a1w/2a1w.pdb1-250.jpg.gz
a1/2a1w/2a1w.pdb1-500.jpg.gz
a1/2a1w/2a1w.pdb1-65.jpg.gz
a1/2a1w/2a1w.pdb1-80.jpg.gz
a1/2a1w/2a1w.pdb1.gz
a1/2a1w/2a1w.pdb2-250.jpg.gz
a1/2a1w/2a1w.pdb2-500.jpg.gz
a1/2a1w/2a1w.pdb2-65.jpg.gz
a1/2a1w/2a1w.pdb2-80.jpg.gz
a1/2a1w/2a1w.pdb2.gz
a1/2a1w/2a1w.pdb3-250.jpg.gz
a1/2a1w/2a1w.pdb3-500.jpg.gz
a1/2a1w/2a1w.pdb3-65.jpg.gz
a1/2a1w/2a1w.pdb3-80.jpg.gz
a1/2a1w/2a1w.pdb3.gz
a1/2a1w/2a1w.pdb4-250.jpg.gz
a1/2a1w/2a1w.pdb4-500.jpg.gz
a1/2a1w/2a1w.pdb4-65.jpg.gz
a1/2a1w/2a1w.pdb4-80.jpg.gz
a1/2a1w/2a1w.pdb4.gz
a1/2a1w/2a1w.pdb5-250.jpg.gz
a1/2a1w/2a1w.pdb5-500.jpg.gz
a1/2a1w/2a1w.pdb5-65.jpg.gz
a1/2a1w/2a1w.pdb5-80.jpg.gz
a1/2a1w/2a1w.pdb5.gz
a1/2a1w/2a1w.pdb6-250.jpg.gz
a1/2a1w/2a1w.pdb6-500.jpg.gz
a1/2a1w/2a1w.pdb6-65.jpg.gz
a1/2a1w/2a1w.pdb6-80.jpg.gz
a1/2a1w/2a1w.pdb6.gz
a1/2a1w/2a1w.stride.gz
a1/2a1w/2a1w.xml.gz
a1/2a1w/pdb2a1w.ent.gz
a1/2a1w/r2a1wsf.ent.gz
a1/2a1x/2a1x-deriv.cif.gz
a1/2a1x/2a1x-extatom.xml.gz
a1/2a1x/2a1x-noatom.xml.gz
</pre>
 * 
 * 
 * @author Andreas Prlic
 *
 *
 *@ since3.2
 */
public class SandboxStyleStructureProvider implements StructureProvider {
	FileParsingParameters params ;
	
	String path;
	public static final String fileSeparator = System.getProperty("file.separator");
	
	public SandboxStyleStructureProvider() {
		params = new FileParsingParameters();
		
		
		String defaultPath = System.getProperty(AbstractUserArgumentProcessor.PDB_DIR);
		if ( defaultPath == null) {
			String property = "java.io.tmpdir";
			defaultPath = System.getProperty(property);
		}
		
		setPath(defaultPath);
	}
	
	/** directory where to find PDB files */
	public void setPath(String p){
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,p);
		
		path = p ;
		
		if ( ! (path.endsWith(fileSeparator) ) )
			path = path + fileSeparator;
				
	}
	
	
	public Structure getStructureById(String pdbId) throws IOException,StructureException {

		
		if (pdbId == null || pdbId.length()< 4)
			throw new StructureException("This does not look like a valid PDB ID! (" + pdbId + ")");
		
		pdbId = pdbId.toLowerCase();
		
		String middle = pdbId.substring(1,3).toLowerCase();
		
		File f = new File(path + fileSeparator + middle + fileSeparator + pdbId  + fileSeparator + "pdb" + pdbId + ".ent.gz");
		
		if (! f.exists()){
			
		}
		
		
		InputStreamProvider isp = new InputStreamProvider();
		
		try {
			InputStream inputStream = isp.getInputStream(f);
			PDBFileParser pdbpars = new PDBFileParser();
			pdbpars.setFileParsingParameters(params);

			Structure struc = pdbpars.parsePDBFile(inputStream) ;
			return struc ;
			
			
		} catch (Exception e){
			e.printStackTrace();
			// something is wrong with the file!
			// it probably should be downloaded again...
			// TODO: add auto-download funcitonality...
		}
		
		return null;
		
		
	}

	public Structure getBiologicalUnit(String pdbId) throws StructureException,
			IOException {
		// TODO Auto-generated method stub
		return null;
	}

	public void setFileParsingParameters(FileParsingParameters params) {
		this.params = params;
	}

	public FileParsingParameters getFileParsingParameters() {
		return params;
	}

	/** Returns a list of all PDB IDs that are available in this installation
	 * 
	 * @return a list of PDB IDs
	 */
	public List<String> getAllPDBIDs() throws IOException{
		
		File f = new File(path);
		if ( ! f.isDirectory())
			throw new IOException("Path " + path + " is not a directory!");
		
		String[] dirName = f.list();
		
		List<String>pdbIds = new ArrayList<String>();
		for (String dir : dirName) {
			File d2= new File(f,dir);
			if ( ! d2.isDirectory())
				continue;
			
			String[] pdbDirs = d2.list();
			for (String pdbId : pdbDirs) {
				if ( ! pdbIds.contains(pdbId))
					pdbIds.add(pdbId);
				
			}
		}
		
		return pdbIds;
	}

}
