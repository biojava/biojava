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
 * created at Sep 7, 2007
 */
package org.biojava.bio.structure.server;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import java.util.logging.Logger;

import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;



/** a PDB installation that reads all files from one directory in the file system.
 * and keeps several files with the information parsed from the PDB headers for
 * easier data access.
 * <br/>
 * In particular these files are:
 * <ul>
 * <li>pdbinfo.txt - contains resolution, experimetnal technique, deposition and 
 *     modification date, biological description of the structure</li>
 * <li>chaininfo.txt - contains information about every chain, like the compound information from the header, the nr of amino acids in chain,
 *     The alignment to the SEQRES residues </li>
 * </ul>
 *     
 * In order to create these index files for this PDB installation, please see the PrepareIndexFile class.
 * @see PrepareIndexFile
 * 
 * @author Andreas Prlic
 * @deprecated
 */
public class FlatFileInstallation 
implements PDBInstallation
{

    public static final Logger logger = Logger.getLogger("org.biojava.bio.structure");

   
    private File filePath;

    private List<PDBFilter> filters;

    public static final String DEFAULT_INDEX_FILE = "pdbinfo.txt";
    public static final String DEFAUL_CHAIN_FILE  = "chaininfo.txt";

    private File indexFile;
    private File chainInfoFile;

    List<PDBHeader> filteredPDBs  ;
    PDBFileReader reader;
    Iterator<PDBHeader> filterIterator;

    boolean filtersApplied ;
    
    public static SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd", Locale.US);
    
    
    /** create a new FlatFile Installation and point it to the directory that contains all PDB files
     * 
     * @param filePath
     */
    public FlatFileInstallation(File filePath){

        if (! filePath.isDirectory()){
            throw new IllegalArgumentException("the provided path does not point to a directory!");
        }


        filters = new ArrayList<PDBFilter>();
        this.filePath = filePath;
        //default is in the same directory as filePath...
        indexFile = new File(filePath + File.separator +DEFAULT_INDEX_FILE);
        chainInfoFile = new File(filePath + File.separator + DEFAUL_CHAIN_FILE);
        filtersApplied = false;
        filteredPDBs = new ArrayList<PDBHeader>();
        filterIterator = null;
        reader = new PDBFileReader();
        reader.setPath(filePath.toString());

    }



    public File getPDBInfoFile() {
        return indexFile;
    }
    public void setPDBInfoFile(File indexFile) {
        this.indexFile = indexFile;
    }
    public File getChainInfoFile() {
        return chainInfoFile;
    }
    public void setChainInfoFile(File chainInfoFile) {
        this.chainInfoFile = chainInfoFile;
    }
    public File getFilePath() {
        return filePath;
    }

    public void setFilePath(File filePath) {
        this.filePath = filePath;
        reader.setPath(filePath.toString());
    }

    public void addPDBFilter(PDBFilter filter) {
        filters.add(filter);
    }

    public void clearFilters() {
        filters.clear();

    }

    /*private void resetIterator() {
        if ( ! filtersApplied )
            applyFilters();
        filterIterator = filteredPDBs.iterator();       
    }*/

    public void applyFilters(){

        filteredPDBs.clear();

        try {

            InputStream in = new FileInputStream(indexFile);
            BufferedReader	indexReader = new BufferedReader(new InputStreamReader(in));

            String line = null;
            line = indexReader.readLine();
            while ( line != null  ) {
                if (line.startsWith("//")) {
                    line = indexReader.readLine();
                    continue;
                }
                               

                PDBHeader header = getPDBHeaderFromLine(line);
                // if no filters have been configured, add all PDBs
                if ( filters.size() == 0){
                    filteredPDBs.add(header);
                } else {
                    throw new UnsupportedOperationException("can not sort PDBHeaders yet");
                }


                line = indexReader.readLine();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }


        filterIterator = filteredPDBs.iterator();		
        filtersApplied = true;
    }

    private String nullCheck(String str){
        if ( str.equals("null")){
            str = null;
        }
        return str;
    }

    private PDBHeader getPDBHeaderFromLine(String line){
        PDBHeader header = new PDBHeader();

        String[] spl = line.split("\t");
        if ( spl.length != 9){
            System.err.println("nr tabs length does not match expected length 9 (is "+spl.length+")! " + line);
        }
//      pdbId\tnrCAAtoms\ttechnique\tresolution\tdepDate\tmodDate\ttitle\tclassification\ttime");
        header.setIdCode(spl[0]);
        
        header.setTechnique(nullCheck(spl[2]));
        try {
            float resolution = Float.parseFloat(spl[3]);
            header.setResolution(resolution);
        } catch (NumberFormatException ex){}

        try {
            
            Date dep = dateFormat.parse(spl[4]);
            header.setDepDate(dep);
            Date mod = dateFormat.parse(spl[5]);
            header.setModDate(mod);
            
        } catch (ParseException e){
            e.printStackTrace();
        }
        
        
        header.setTitle(nullCheck(spl[6]));
        header.setClassification(nullCheck(spl[7]));

        return header;
    }

    public List<PDBHeader> getAll() {

        if ( ! filtersApplied)
            applyFilters();

        throw new UnsupportedOperationException("can not get all PDBHeaders yet");

    }

    public PDBHeader getPDBHeader(String pdbId){

        Iterator<PDBHeader> iter = filteredPDBs.iterator();

        while (iter.hasNext()){
            PDBHeader header = iter.next();
            String id = header.getIdCode();
            if ( pdbId.equals(id)){
                return header;
            }    
        }

        return null;
    }

    public Structure getStructure(String pdbId) {
        Structure s = null;
        try {
            s = reader.getStructureById(pdbId);
        } catch (IOException e){
            e.printStackTrace();
        }
        return s;
    }

    public boolean hasNext() {
        if ( ! filtersApplied)
            applyFilters();

        return filterIterator.hasNext(); 
    }


    public Structure next() {


        if ( ! filtersApplied )
            applyFilters();

        PDBHeader header = filterIterator.next();
        String pdbId = header.getIdCode();

        Structure s = null;
        try {
            System.out.println("flatfileinstallation : next:" + pdbId + " " + getFilePath());
            s = reader.getStructureById(pdbId);
        } catch (IOException e){
            e.printStackTrace();
        }
        return s;


    }



}
